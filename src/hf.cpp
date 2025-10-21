#include "hf.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

namespace hf {

static int electron_count(const Molecule& mol) {
  int Zsum = 0;
  for (const auto& a : mol.atoms()) Zsum += a.Z;
  return Zsum - mol.net_charge();
}

HF::HF(std::shared_ptr<Molecule> mol, double conv_thresh, size_t max_iter)
  : mol_(std::move(mol)), conv_thresh_(conv_thresh), max_iter_(max_iter) {
  const size_t n = nbf();
  F_.setZero(n, n);
  C_.setZero(n, n);
  R_.setZero(n, n);
  J_.setZero(n, n);
  K_.setZero(n, n);
  G_.setZero(n, n);
}

// Build density for RHF: R = 2 sum_i occ C(:,i) C(:,i)^T
static void build_density_RHF(const Eigen::MatrixXd& C, int nelec, Eigen::MatrixXd& R) {
  const auto nbf = C.rows();
  R.setZero(nbf, nbf);
  int nocc = nelec / 2; // closed-shell
  Eigen::MatrixXd Cocc = C.leftCols(nocc);
  R = 2.0 * (Cocc * Cocc.transpose());
}

// Build J and K from ERI and density
static inline size_t eri_index(size_t mu, size_t nu, size_t lam, size_t sig, size_t n) {
  return ((mu * n + nu) * n + lam) * n + sig;
}
static void build_JK(const std::vector<double>& eri, const Eigen::MatrixXd& R, size_t nbf,
                     Eigen::MatrixXd& J, Eigen::MatrixXd& K) {
  J.setZero(nbf, nbf);
  K.setZero(nbf, nbf);
  for (size_t mu = 0; mu < nbf; ++mu) {
    for (size_t nu = 0; nu < nbf; ++nu) {
      double Jmunu = 0.0;
      double Kmunu = 0.0;
      for (size_t lam = 0; lam < nbf; ++lam) {
        for (size_t sig = 0; sig < nbf; ++sig) {
          double P = R(lam, sig);
          Jmunu += P * eri[eri_index(mu, nu, lam, sig, nbf)];
          Kmunu += P * eri[eri_index(mu, lam, nu, sig, nbf)];
        }
      }
      J(mu, nu) = Jmunu;
      K(mu, nu) = Kmunu;
    }
  }
}

// Energy: E = sum_{μν} R_{μν} (h_{μν} + F_{μν})/2 = tr( R (h + F) ) / 2
double HF::compute_energy() const {
  auto h = const_cast<HF*>(this)->get_hcore();
  Eigen::MatrixXd sumHF = h + F_;
  return 0.5 * (R_.cwiseProduct(sumHF)).sum();
}

bool HF::check_convergence() const {
  // simple: density RMS change and energy change threshold
  if (iter_ == 0) return false;
  double dE = std::abs(energy_ - prev_energy_);
  double dR_rms = (R_ - prev_R_).norm() / std::sqrt(static_cast<double>(R_.size()));
  return (dE < conv_thresh_) && (dR_rms < conv_thresh_);
}

bool HF::iterate() {
  const size_t n = nbf();
  // cache previous
  prev_R_ = R_;
  prev_energy_ = energy_;

  // Build Fock: F = h + G, with G = 2J - K (RHF)
  auto h = get_hcore();
  auto eri = get_eri();
  build_JK(eri, R_, n, J_, K_);
  G_ = 2.0 * J_ - K_;
  F_ = h + G_;

  // Solve generalized eigenproblem F C = S C e
  // Build Eigen matrices
  Eigen::MatrixXd F = F_;
  Eigen::MatrixXd S = mol_->compute_overlap();

  // Symmetric orthogonalization: S = U s U^T, X = U s^{-1/2} U^T
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
  Eigen::MatrixXd U = es.eigenvectors();
  Eigen::VectorXd s = es.eigenvalues();
  Eigen::VectorXd s_inv_sqrt = s.array().inverse().sqrt();
  Eigen::MatrixXd Xinvsqrt = U * s_inv_sqrt.asDiagonal() * U.transpose();

  // Transform F to orthonormal basis and diagonalize
  Eigen::MatrixXd Fp = Xinvsqrt.transpose() * F * Xinvsqrt;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esF(Fp);
  Eigen::MatrixXd Cp = esF.eigenvectors();
  Eigen::MatrixXd C = Xinvsqrt * Cp; // back-transform to AO basis

  // Update C_
  C_ = C;

  // Update density from occupied orbitals
  int nelec = electron_count(*mol_);
  build_density_RHF(C_, nelec, R_);

  // Update energy
  energy_ = compute_energy();
  ++iter_;
  converged_ = check_convergence();
  return converged_;
}

void HF::run() {
  // Initialize: start from core Hamiltonian diagonalization density
  const size_t n = nbf();
  // Initial Fock = hcore
  auto h = get_hcore();
  F_ = h;

  // Overlap
  Eigen::MatrixXd F = F_;
  Eigen::MatrixXd S = mol_->compute_overlap();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
  Eigen::MatrixXd U = es.eigenvectors();
  Eigen::VectorXd s = es.eigenvalues();
  Eigen::VectorXd s_inv_sqrt = s.array().inverse().sqrt();
  Eigen::MatrixXd Xinvsqrt = U * s_inv_sqrt.asDiagonal() * U.transpose();
  Eigen::MatrixXd Fp = Xinvsqrt.transpose() * F * Xinvsqrt;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esF(Fp);
  Eigen::MatrixXd Cp = esF.eigenvectors();
  Eigen::MatrixXd C = Xinvsqrt * Cp;
  C_ = C;

  int nelec = electron_count(*mol_);
  build_density_RHF(C_, nelec, R_);
  energy_ = compute_energy();
  prev_energy_ = energy_;
  prev_R_ = R_;

  for (iter_ = 0; iter_ < max_iter_; ++iter_) {
    if (iterate()) break;
  }
}

} // namespace hf

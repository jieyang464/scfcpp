// Hartree-Fock state and SCF loop for a molecule
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cstddef>
#include <memory>

namespace hf {

// Forward declaration to avoid including molecule.h here
class Molecule;

class HF {
public:
  HF(std::shared_ptr<Molecule> mol, double conv_thresh = 1e-8, size_t max_iter = 50);
  ~HF(); // out-of-line dtor so Molecule is complete at destruction

  // Run SCF loop
  void run();
  // Perform one SCF iteration
  bool iterate();
  // Compute energy
  double compute_energy() const;
  // Check convergence
  bool check_convergence() const;

  // Accessors
  const Eigen::MatrixXd& get_F() const { return F_; }
  const Eigen::MatrixXd& get_C() const { return C_; }
  const Eigen::MatrixXd& get_R() const { return R_; }
  const Eigen::MatrixXd& get_J() const { return J_; }
  const Eigen::MatrixXd& get_K() const { return K_; }
  const Eigen::MatrixXd& get_G() const { return G_; }
  double energy() const { return energy_; }
  bool converged() const { return converged_; }
  size_t nbf() const;
  const Molecule& molecule() const;
  int nelec_alpha() const { return nelec_alpha_; }
  int nelec_beta()  const { return nelec_beta_;  }

private:
  std::shared_ptr<Molecule> mol_;
  // RHF matrices
  Eigen::MatrixXd F_, C_, R_, J_, K_, G_;
  double energy_ = 0.0;
  double prev_energy_ = 0.0; // used for convergence check
  bool converged_ = false;
  double conv_thresh_ = 1e-8;
  size_t max_iter_ = 50;
  size_t iter_ = 0;
  int nelec_alpha_ = 0;
  int nelec_beta_  = 0;
  Eigen::MatrixXd prev_R_; // used for convergence check

  // Helper: get hcore and ERI from molecule
  Eigen::MatrixXd get_hcore() const;
  std::vector<double> get_eri() const;
};

} // namespace hf
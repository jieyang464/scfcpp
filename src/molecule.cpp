// Implementation for Molecule using libint2

#include "molecule.h"
#include <libint2.hpp>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>

namespace hf {
// RAII guard to ensure libint2 is initialized before use and finalized at program end
struct LibintGuard {
  LibintGuard() { libint2::initialize(); }
  ~LibintGuard() { libint2::finalize(); }
};

// Define proxy struct names to avoid exposing libint types in header
struct libint2_Atom_proxy { int dummy; };
struct libint2_BasisSet_proxy { int dummy; };

// forward declaration of helper
static libint2::BasisSet choose_basis(const Molecule& mol, const std::vector<libint2::Atom>& ai);

std::vector<libint2_Atom_proxy> Molecule::to_libint_atoms() const {
  std::vector<libint2_Atom_proxy> proxy; proxy.reserve(atoms_.size());
  // We will build actual libint2::Atom on the fly when needed
  for (size_t i = 0; i < atoms_.size(); ++i) proxy.push_back({0});
  return proxy;
}

static std::vector<libint2::Atom> build_libint_atoms(const std::vector<Atom>& atoms) {
  std::vector<libint2::Atom> ai; ai.reserve(atoms.size());
  for (const auto& a : atoms) {
    libint2::Atom la; la.atomic_number = a.Z; la.x = a.x; la.y = a.y; la.z = a.z; ai.push_back(la);
  }
  return ai;
}

libint2_BasisSet_proxy Molecule::make_basis() const {
  return {0};
}

size_t Molecule::nbf() const {
  static LibintGuard guard; // initialize once per process
  auto ai = build_libint_atoms(atoms_);
  libint2::BasisSet basis = choose_basis(*this, ai);
  return basis.nbf();
}

// Helper to choose basis: if basis_file_path_ empty, check ./basis/<basis_name>.g94 relative path first
static libint2::BasisSet choose_basis(const Molecule& mol, const std::vector<libint2::Atom>& ai) {
  // If user provided an explicit file path, set LIBINT_DATA_PATH to its directory and use basename
  if (!mol.basis_file_path().empty()) {
    std::string path = mol.basis_file_path();
    auto pos = path.find_last_of('/');
    std::string dir = (pos==std::string::npos) ? std::string(".") : path.substr(0,pos);
    std::string file = (pos==std::string::npos) ? path : path.substr(pos+1);
    std::string name = file;
    // strip .g94 if present
    if (name.size() > 4 && name.substr(name.size()-4)==".g94") name = name.substr(0, name.size()-4);
    setenv("LIBINT_DATA_PATH", dir.c_str(), 1);
    return libint2::BasisSet(name, ai, /*throw_if_no_match=*/true);
  }

  // try relative path ./basis/<name>.g94
  std::string rel = std::string("basis/") + mol.basis_name() + ".g94";
  if (std::ifstream(rel)) {
    setenv("LIBINT_DATA_PATH", "./basis", 1);
    return libint2::BasisSet(mol.basis_name(), ai, /*throw_if_no_match=*/true);
  }

  // fallback to built-in (system) data path
  return libint2::BasisSet(mol.basis_name(), ai, /*throw_if_no_match=*/false);
}

Eigen::MatrixXd Molecule::compute_overlap() const {
  static LibintGuard guard;
  return compute_one_electron(static_cast<int>(libint2::Operator::overlap));
}

Eigen::MatrixXd Molecule::compute_kinetic() const {
  static LibintGuard guard;
  return compute_one_electron(static_cast<int>(libint2::Operator::kinetic));
}

Eigen::MatrixXd Molecule::compute_nuclear_attraction() const {
  static LibintGuard guard;
  auto ai = build_libint_atoms(atoms_);
  libint2::BasisSet basis = choose_basis(*this, ai);
  const size_t nbasis = basis.nbf();
  Eigen::MatrixXd V(nbasis, nbasis);
  V.setZero();

  libint2::Engine V_engine(libint2::Operator::nuclear, basis.max_nprim(), basis.max_l());
  V_engine.set_params(libint2::make_point_charges(ai));
  const auto& Vbufs = V_engine.results();

  size_t row_offset = 0;
  for (size_t s1 = 0; s1 != basis.size(); ++s1) {
    const auto& sh1 = basis[s1];
    const size_t n1 = sh1.size();
    size_t col_offset = 0;
    for (size_t s2 = 0; s2 != basis.size(); ++s2) {
      const auto& sh2 = basis[s2];
      const size_t n2 = sh2.size();
      V_engine.compute(sh1, sh2);
      const auto* buf = Vbufs[0];
      if (buf == nullptr) { col_offset += n2; continue; }
      for (size_t f1 = 0; f1 != n1; ++f1) {
        for (size_t f2 = 0; f2 != n2; ++f2) {
          V(row_offset + f1, col_offset + f2) = buf[f1 * n2 + f2];
        }
      }
      col_offset += n2;
    }
    row_offset += n1;
  }
  return V;
}

double Molecule::compute_eri_abs_sum() const {
  static LibintGuard guard;
  auto ai = build_libint_atoms(atoms_);
  libint2::BasisSet basis = choose_basis(*this, ai);
  libint2::Engine eri_engine(libint2::Operator::coulomb, basis.max_nprim(), basis.max_l());
  const auto& ERIbufs = eri_engine.results();

  double eri_sum = 0.0;
  for (size_t s1 = 0; s1 != basis.size(); ++s1) {
    const auto& sh1 = basis[s1];
    for (size_t s2 = 0; s2 != basis.size(); ++s2) {
      const auto& sh2 = basis[s2];
      for (size_t s3 = 0; s3 != basis.size(); ++s3) {
        const auto& sh3 = basis[s3];
        for (size_t s4 = 0; s4 != basis.size(); ++s4) {
          const auto& sh4 = basis[s4];
          eri_engine.compute(sh1, sh2, sh3, sh4);
          const auto* buf = ERIbufs[0];
          if (!buf) continue;
          const size_t n1 = sh1.size(), n2 = sh2.size();
          const size_t n3 = sh3.size(), n4 = sh4.size();
          const size_t n1234 = n1 * n2 * n3 * n4;
          for (size_t i = 0; i < n1234; ++i) eri_sum += std::abs(buf[i]);
        }
      }
    }
  }
  return eri_sum;
}

Eigen::MatrixXd Molecule::compute_one_electron(int op_int) const {
  static LibintGuard guard;
  libint2::Operator op = static_cast<libint2::Operator>(op_int);
  auto ai = build_libint_atoms(atoms_);
  libint2::BasisSet basis = choose_basis(*this, ai);
  const size_t nbasis = basis.nbf();
  Eigen::MatrixXd M(nbasis, nbasis);
  M.setZero();

  libint2::Engine engine(op, basis.max_nprim(), basis.max_l());
  const auto& bufs = engine.results();

  size_t row_offset = 0;
  for (size_t s1 = 0; s1 != basis.size(); ++s1) {
    const auto& sh1 = basis[s1];
    const size_t n1 = sh1.size();
    size_t col_offset = 0;
    for (size_t s2 = 0; s2 != basis.size(); ++s2) {
      const auto& sh2 = basis[s2];
      const size_t n2 = sh2.size();
      engine.compute(sh1, sh2);
      const auto* buf = bufs[0];
      if (buf == nullptr) { col_offset += n2; continue; }
      for (size_t f1 = 0; f1 != n1; ++f1) {
        for (size_t f2 = 0; f2 != n2; ++f2) {
          M(row_offset + f1, col_offset + f2) = buf[f1 * n2 + f2];
        }
      }
      col_offset += n2;
    }
    row_offset += n1;
  }
  return M;
}

std::vector<double> Molecule::compute_eri_tensor() const {
  static LibintGuard guard;
  auto ai = build_libint_atoms(atoms_);
  libint2::BasisSet basis = choose_basis(*this, ai);
  const size_t nbf = basis.nbf();
  std::vector<double> eri(nbf * nbf * nbf * nbf, 0.0);
  auto eri_index = [nbf](size_t mu, size_t nu, size_t lam, size_t sig) {
    return ((mu * nbf + nu) * nbf + lam) * nbf + sig;
  };

  libint2::Engine eri_engine(libint2::Operator::coulomb, basis.max_nprim(), basis.max_l());
  const auto& ERIbufs = eri_engine.results();
  for (size_t s1 = 0; s1 != basis.size(); ++s1) {
    const auto& sh1 = basis[s1];
    const size_t n1 = sh1.size();
    for (size_t s2 = 0; s2 != basis.size(); ++s2) {
      const auto& sh2 = basis[s2];
      const size_t n2 = sh2.size();
      for (size_t s3 = 0; s3 != basis.size(); ++s3) {
        const auto& sh3 = basis[s3];
        const size_t n3 = sh3.size();
        for (size_t s4 = 0; s4 != basis.size(); ++s4) {
          const auto& sh4 = basis[s4];
          const size_t n4 = sh4.size();
          eri_engine.compute(sh1, sh2, sh3, sh4);
          const auto* buf = ERIbufs[0];
          if (!buf) continue;
          for (size_t f1 = 0; f1 != n1; ++f1) {
            for (size_t f2 = 0; f2 != n2; ++f2) {
              for (size_t f3 = 0; f3 != n3; ++f3) {
                for (size_t f4 = 0; f4 != n4; ++f4) {
                  size_t mu = /* row offset */ 0; // compute offsets
                }
              }
            }
          }
          // compute offsets to place buf into eri array
          // calculate basis function offsets
          static size_t base_row_offsets[1]; // dummy
          // We'll compute function offsets by traversing shells again
          // Build offsets
          std::vector<size_t> shell_offsets(basis.size()+1, 0);
          for (size_t s=0; s < basis.size(); ++s) shell_offsets[s+1] = shell_offsets[s] + basis[s].size();
          size_t off1 = shell_offsets[s1];
          size_t off2 = shell_offsets[s2];
          size_t off3 = shell_offsets[s3];
          size_t off4 = shell_offsets[s4];
          for (size_t f1 = 0; f1 != n1; ++f1)
            for (size_t f2 = 0; f2 != n2; ++f2)
              for (size_t f3 = 0; f3 != n3; ++f3)
                for (size_t f4 = 0; f4 != n4; ++f4) {
                  size_t mu = off1 + f1;
                  size_t nu = off2 + f2;
                  size_t lam = off3 + f3;
                  size_t sig = off4 + f4;
                  size_t idx = eri_index(mu, nu, lam, sig);
                  size_t buf_idx = ((f1 * n2 + f2) * n3 + f3) * n4 + f4; // keep local layout for buf
                  eri[idx] = buf[buf_idx];
                }
        }
      }
    }
  }
  return eri;
}

Eigen::MatrixXd Molecule::compute_hcore() const {
  static LibintGuard guard;
  auto T = compute_kinetic();
  auto V = compute_nuclear_attraction();
  if (T.rows() != V.rows() || T.cols() != V.cols()) throw std::runtime_error("T and V dimension mismatch");
  Eigen::MatrixXd H = T + V;
  return H;
}

Eigen::MatrixXd Molecule::evaluate_ao(const Eigen::MatrixXd& points) const {
  // Placeholder: Implement AO evaluation using libint or a custom evaluator.
  // Contract: returns (n_points x nbf) matrix, columns ordered by atom order following BasisSet.
  // Current behavior: throw until implemented to avoid silent misuse.
  throw std::runtime_error("evaluate_ao is not implemented yet. Provide a libint-based AO evaluator.");
}

} // namespace hf

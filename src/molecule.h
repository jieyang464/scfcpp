  // (AO grid evaluation API removed — implement later using a dedicated AO evaluator.)
// Simple Molecule/Atom definitions and libint2-based integral computations
// Assumptions:
// - A single uniform basis set is used for all atoms in the molecule (by name and optional file path).
// - Outputs for 1e operators (S, T, V) are returned as row-major vectors of size nbf*nbf.
// - ERI computation provides a simple aggregated measure or can be extended to build a full 4-index tensor.
// - Molecule internally manages libint2 initialization/finalization; no external calls needed.

#pragma once

#include <vector>
#include <string>
#include <Eigen/Dense>

namespace hf {

struct Atom {
  int Z;           // atomic number
  double x, y, z;  // Cartesian coordinates in Angstrom
  // Optional: per-atom basis label (not used in BasisSet construction here; kept for future extension)
  std::string basis_label;
};

class Molecule {
public:
  Molecule() = default;

  Molecule(std::vector<Atom> atoms, int net_charge,
           std::string basis_name,
           std::string basis_file_path = "")
      : atoms_(std::move(atoms)), net_charge_(net_charge),
        basis_name_(std::move(basis_name)), basis_file_path_(std::move(basis_file_path)) {}

  // Accessors
  const std::vector<Atom>& atoms() const { return atoms_; }
  int net_charge() const { return net_charge_; }
  const std::string& basis_name() const { return basis_name_; }
  const std::string& basis_file_path() const { return basis_file_path_; }

  // Build libint2 atoms vector
  std::vector<struct libint2_Atom_proxy> to_libint_atoms() const; // implemented in .cpp via libint2::Atom

  // Construct BasisSet; if basis_file_path_ is provided, use it, otherwise rely on built-in data
  // Construct BasisSet; implemented in .cpp
  struct libint2_BasisSet_proxy make_basis() const; // implemented in .cpp via libint2::BasisSet

  // Number of basis functions for the molecule
  size_t nbf() const;

  // Compute overlap matrix (nbf x nbf)
  Eigen::MatrixXd compute_overlap() const;

  // Compute kinetic energy matrix (nbf x nbf)
  Eigen::MatrixXd compute_kinetic() const;

  // Compute nuclear attraction matrix (sum of -Z/|r-R_A| over all nuclei)
  Eigen::MatrixXd compute_nuclear_attraction() const;

  // Compute ERIs: returns a simple accumulated absolute sum for sanity check.
  // For production, consider building and storing a 4-index tensor (μν|λσ) with permutational symmetry.
  double compute_eri_abs_sum() const;

  // Compute full ERI tensor in chemists' notation (mu,nu,lam,sig) as a flat vector of size nbf^4
  // Indexing: idx = ((mu * nbf + nu) * nbf + lam) * nbf + sig
  std::vector<double> compute_eri_tensor() const;

  // Compute core Hamiltonian H = T + V (nbf x nbf)
  Eigen::MatrixXd compute_hcore() const;

  // ---------------------------------------------------------------------------
  // SCF related types and API
  // The Molecule class exposes a simple SCF entry point that runs an RHF solver
  // for this geometry/basis. The solver implementation lives in the .cpp and is
  // internal to the hf module; the method returns a small result container and
  // also stores the last-converged 1-RDM / canonical MOs in the Molecule instance
  // for convenient access by geometry optimization and downstream code.
  struct HFResult {
    double energy = 0.0;                 // total electronic energy (without nuclear repulsion)
    Eigen::MatrixXd density;             // AO 1-RDM (nbf x nbf)
    Eigen::MatrixXd C;                   // MO coefficients (nbf x nbf)
    Eigen::VectorXd eps;                 // orbital energies (nbf)
    bool converged = false;
    int niter = 0;
  };

  // Run a simple RHF SCF for this molecule. The implementation uses the
  // molecule's basis/geometry and returns HFResult; it will also populate
  // the molecule's internal cached density/C/eps so callers can access them
  // without threading HFResult around. Default options are conservative; add
  // optional parameters in the .cpp as needed.
  HFResult run_scf(double conv_tol = 1e-8, int max_iter = 100);

  // Accessors for the last computed SCF quantities. They will be empty if no
  // SCF has been run yet for this Molecule instance.
  bool has_last_density() const { return has_last_density_; }
  const Eigen::MatrixXd& last_density() const { return last_density_; }
  const Eigen::MatrixXd& last_C() const { return last_C_; }
  const Eigen::VectorXd& last_eps() const { return last_eps_; }
  double last_scf_energy() const { return last_energy_; }
  int last_scf_iterations() const { return last_scf_iterations_; }
  bool last_scf_converged() const { return last_scf_converged_; }
  // ---------------------------------------------------------------------------

  // Evaluate AO basis functions on a set of 3D points.
  // Returns a matrix with shape (n_points x nbf), where columns correspond to AO basis functions
  // ordered by atom order in the molecule, expanding shells and Cartesian/GTO primitives consistently
  // with the libint BasisSet ordering.
  // points: each row is (x,y,z) in Angstrom.
  Eigen::MatrixXd evaluate_ao(const Eigen::MatrixXd& points) const;

private:
  // Helper to compute generic one-electron operator
  Eigen::MatrixXd compute_one_electron(int op) const; // implemented in .cpp using libint2::Operator

  std::vector<Atom> atoms_;
  int net_charge_ = 0;
  std::string basis_name_ = "sto-3g";
  std::string basis_file_path_;
};

} // namespace hf

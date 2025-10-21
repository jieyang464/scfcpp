#include "molecule.h"
#include <iostream>
#include <iomanip>

int main() {

  // HeH+: two atoms separated by 1.4632 a.u.
  // Place He at origin and H along z at 1.4632
  hf::Atom He{2, 0.0, 0.0, 0.0};
  hf::Atom H{1, 0.0, 0.0, 1.4632};
  hf::Molecule mol({He, H}, +1, "sto-3g", ""); // empty path -> uses relative basis/built-in

  try {
  Eigen::MatrixXd hcore = mol.compute_hcore();
    auto eri = mol.compute_eri_tensor();
    size_t nbf = mol.nbf();
    std::cout << "nbf=" << nbf << "\n";

    std::cout << "Hcore (row-major):\n";
    for (size_t i = 0; i < nbf; ++i) {
      for (size_t j = 0; j < nbf; ++j) std::cout << std::setw(15) << hcore(i, j);
      std::cout << '\n';
    }

    std::cout << "\nERI tensor (non-zero list)\n";
    // print only unique ERIs where mu<=nu and lam<=sig to reduce output
    for (size_t mu = 0; mu < nbf; ++mu)
      for (size_t nu = 0; nu < nbf; ++nu)
        for (size_t lam = 0; lam < nbf; ++lam)
          for (size_t sig = 0; sig < nbf; ++sig) {
            size_t idx = ((mu * nbf + nu) * nbf + lam) * nbf + sig;
            double v = eri[idx];
            if (std::abs(v) > 1e-12) {
              std::cout << "("<<mu<<","<<nu<<"|"<<lam<<","<<sig<<") = "<< v <<"\n";
            }
          }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
  return 0;
}

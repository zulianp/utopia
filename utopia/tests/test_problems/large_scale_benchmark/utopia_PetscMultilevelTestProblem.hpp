#ifndef UTOPIA_PETSC_MULTILEVEL_TEST_PROBLEM_HPP
#define UTOPIA_PETSC_MULTILEVEL_TEST_PROBLEM_HPP

#include "utopia.hpp"
#include "utopia_TestFunctions.hpp"

#ifdef WITH_PETSC
#include <petsc/private/snesimpl.h> /* For SNES_Solve event */
#include <petscdm.h>
#include <petscdmda.h>
#include <petscmatlab.h>
#include <petscsnes.h>

namespace utopia {

template <typename Matrix, typename Vector, typename ProblemType>
class PetscMultilevelTestProblem final
    : public MultilevelTestProblemBase<Matrix, Vector> {
 public:
  using SizeType = typename utopia::Traits<Vector>::SizeType;
  using Scalar = typename utopia::Traits<Vector>::Scalar;

  PetscMultilevelTestProblem(const SizeType dimension,
                             const SizeType &n_levels = 2,
                             const SizeType &n_coarse = 10,
                             const bool remove_bc = false)
      : MultilevelTestProblemBase<Matrix, Vector>(n_levels, n_coarse,
                                                  remove_bc) {
    std::vector<DM> dms_;
    dms_.resize(n_levels);
    // level_functions_.resize(n_levels);

    if (dimension == 2) {
      DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                   DMDA_STENCIL_STAR, n_coarse, n_coarse, PETSC_DECIDE,
                   PETSC_DECIDE, 1, 1, nullptr, nullptr, &dms_[0]);
    } else if (dimension == 3) {
      DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                   DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, n_coarse, n_coarse,
                   n_coarse, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1,
                   nullptr, nullptr, nullptr, &dms_[0]);
      // DMDASetInterpolationType(dms_[0], DMDA_Q0);
    } else {
      utopia_error(
          "PetscMultilevelTestProblem:: choose valid dimension (2, 3). ");
    }

    DMSetUp(dms_[0]);
    DMDASetUniformCoordinates(dms_[0], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

    for (auto l = 1; l < n_levels; l++) {
      DMRefine(dms_[l - 1], PETSC_COMM_WORLD, &dms_[l]);
    }

    for (auto l = 1; l < n_levels; l++) {
      Mat I;
      Matrix I_u;
      DMCreateInterpolation(dms_[l - 1], dms_[l], &I, nullptr);
      wrap(I, I_u);
      // transfers_.push_back( std::make_shared<IPTransfer<Matrix, Vector> >(
      // std::make_shared<Matrix>(I_u)) );
      this->transfers_[l - 1] = std::make_shared<IPRTransfer<Matrix, Vector> >(
          std::make_shared<Matrix>(I_u));
      MatDestroy(&I);
    }

    for (auto l = 0; l < n_levels; l++) {
      auto fun = std::make_shared<ProblemType>(dms_[l]);
      this->level_functions_[l] = fun;
    }
  }
};
}  // namespace utopia

#endif  // WITH_PETSC

#endif  // UTOPIA_PETSC_MULTILEVEL_TEST_PROBLEM_HPP

#include "utopia_petsc_LinearSolverFactory.hpp"

#include "utopia_AlgebraicMultigrid.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_Traits.hpp"

#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_KSPSolvers.hpp"
#include "utopia_petsc_RowView.hpp"

#include "utopia_FactoryMethod.hpp"
#include "utopia_petsc_BDDLinearSolver.hpp"
#include "utopia_petsc_Factorization.hpp"
#include "utopia_petsc_Factorizations.hpp"
#include "utopia_petsc_GMRES.hpp"

#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

#include "utopia_make_unique.hpp"
#include "utopia_petsc.hpp"

#include "utopia_petsc_KSPSolverMF.hpp"

#include <map>
#include <memory>
#include <string>

namespace utopia {

    void MatrixFreeLinearSolverFactory<PetscVector, PETSC>::register_solvers() {
        Super::register_solvers();
        this->register_solver<KSP_MF<PetscMatrix, PetscVector, PETSC>>(Solver::ksp());
    }

    std::unique_ptr<LinearSolver<PetscMatrix, PetscVector>>
    LinearSolverFactory<PetscMatrix, PetscVector, PETSC>::new_linear_solver(const std::string &tag) {
        auto it = instance().solvers_.find(tag);
        if (it == instance().solvers_.end()) {
            return default_linear_solver();
        } else {
            return it->second->make();
        }
    }

    std::unique_ptr<LinearSolver<PetscMatrix, PetscVector>>
    LinearSolverFactory<PetscMatrix, PetscVector, PETSC>::default_linear_solver() {
        return utopia::make_unique<GMRES<PetscMatrix, PetscVector>>("bjacobi");
    }

    LinearSolverFactory<PetscMatrix, PetscVector, PETSC>
        &LinearSolverFactory<PetscMatrix, PetscVector, PETSC>::instance() {
        static LinearSolverFactory instance_;
        return instance_;
    }

    LinearSolverFactory<PetscMatrix, PetscVector, PETSC>::LinearSolverFactory() { init(); }

    void LinearSolverFactory<PetscMatrix, PetscVector, PETSC>::init() {
        solvers_[Solver::utopia_cg()] =
            utopia::make_unique<LSFactoryMethod<ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE>>>();
        solvers_[Solver::automatic()] = utopia::make_unique<LSFactoryMethod<BiCGStab<PetscMatrix, PetscVector>>>();
        solvers_[Solver::cg()] = utopia::make_unique<LSFactoryMethod<ConjugateGradient<PetscMatrix, PetscVector>>>();
        solvers_[Solver::bicgstab()] = utopia::make_unique<LSFactoryMethod<BiCGStab<PetscMatrix, PetscVector>>>();
        solvers_[Solver::ksp()] = utopia::make_unique<LSFactoryMethod<KSPSolver<PetscMatrix, PetscVector>>>();
        solvers_[Solver::direct()] = utopia::make_unique<LSFactoryMethod<Factorization<PetscMatrix, PetscVector>>>();
        solvers_[Solver::gmres()] = utopia::make_unique<LSFactoryMethod<GMRES<PetscMatrix, PetscVector>>>();
        solvers_["amg"] = utopia::make_unique<LSFactoryMethod<AlgebraicMultigrid<PetscMatrix, PetscVector>>>();
        solvers_["bdd"] = utopia::make_unique<LSFactoryMethod<BDDLinearSolver<PetscMatrix, PetscVector>>>();
    }

}  // namespace utopia

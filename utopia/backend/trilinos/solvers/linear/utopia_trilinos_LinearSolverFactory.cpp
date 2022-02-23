#include "utopia_trilinos_LinearSolverFactory.hpp"

#include "utopia_AlgebraicMultigrid.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_FactoryMethod.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_Traits.hpp"

#include "utopia_make_unique.hpp"
#include "utopia_trilinos.hpp"

#include <map>
#include <memory>
#include <string>

namespace utopia {

    std::unique_ptr<LinearSolver<TpetraMatrix, TpetraVector>>
    LinearSolverFactory<TpetraMatrix, TpetraVector, TRILINOS>::new_linear_solver(const std::string &tag) {
        auto it = instance().solvers_.find(tag);
        if (it == instance().solvers_.end()) {
            return utopia::make_unique<ConjugateGradient<TpetraMatrix, TpetraVector>>();
        } else {
            return it->second->make();
        }
    }

    std::unique_ptr<LinearSolver<TpetraMatrix, TpetraVector>>
    LinearSolverFactory<TpetraMatrix, TpetraVector, TRILINOS>::default_linear_solver() {
        return utopia::make_unique<ConjugateGradient<TpetraMatrix, TpetraVector>>();
    }

    LinearSolverFactory<TpetraMatrix, TpetraVector, TRILINOS>
        &LinearSolverFactory<TpetraMatrix, TpetraVector, TRILINOS>::instance() {
        static LinearSolverFactory instance_;
        return instance_;
    }

    LinearSolverFactory<TpetraMatrix, TpetraVector, TRILINOS>::LinearSolverFactory() { init(); }

    void LinearSolverFactory<TpetraMatrix, TpetraVector, TRILINOS>::init() {
        solvers_[Solver::utopia_cg()] =
            utopia::make_unique<LSFactoryMethod<ConjugateGradient<TpetraMatrix, TpetraVector, HOMEMADE>>>();
        // solvers_[Solver::automatic()] = utopia::make_unique<LSFactoryMethod<BiCGStab<TpetraMatrix, TpetraVector>>>();
        solvers_[Solver::cg()] = utopia::make_unique<LSFactoryMethod<ConjugateGradient<TpetraMatrix, TpetraVector>>>();
        // solvers_[Solver::bicgstab()] = utopia::make_unique<LSFactoryMethod<BiCGStab<TpetraMatrix, TpetraVector>>>();
        // solvers_[Solver::ksp()] = utopia::make_unique<LSFactoryMethod<KSPSolver<TpetraMatrix, TpetraVector>>>();
        // solvers_[Solver::direct()] = utopia::make_unique<LSFactoryMethod<Factorization<TpetraMatrix,
        // TpetraVector>>>(); solvers_["gmres"] = utopia::make_unique<LSFactoryMethod<GMRES<TpetraMatrix,
        // TpetraVector>>>(); solvers_["amg"] = utopia::make_unique<LSFactoryMethod<AlgebraicMultigrid<TpetraMatrix,
        // TpetraVector>>>();
    }

}  // namespace utopia

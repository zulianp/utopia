#ifndef UTOPIA_BLAS_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_BLAS_LINEAR_SOLVER_FACTORY_HPP

#include "utopia_Base.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_FactoryMethod.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_MatrixFreeLinearSolverFactory.hpp"
#include "utopia_Traits.hpp"
#include "utopia_make_unique.hpp"

#ifdef UTOPIA_ENABLE_LAPACK
#include "utopia_Lapack.hpp"
#endif  // UTOPIA_ENABLE_LAPACK

#ifdef UTOPIA_ENABLE_UMFPACK
#include "utopia_UmfpackLU.hpp"
#endif  // UTOPIA_ENABLE_UMFPACK

#include <map>
#include <memory>
#include <string>

namespace utopia {

    template <typename Vector>
    class MatrixFreeLinearSolverFactory<Vector, BLAS> : public BasicMatrixFreeLinearSolverFactory<Vector> {};

    template <typename Matrix, typename Vector>
    class LinearSolverFactory<Matrix, Vector, BLAS> {
    public:
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolverT;
        using LinearSolverPtr = std::unique_ptr<LinearSolverT>;
        using FactoryMethodT = utopia::IFactoryMethod<LinearSolverT>;

        template <class Alg>
        using LSFactoryMethod = FactoryMethod<LinearSolverT, Alg>;
        std::map<std::string, std::shared_ptr<FactoryMethodT>> solvers_;

        inline static LinearSolverPtr new_linear_solver(const std::string &tag) {
            auto it = instance().solvers_.find(tag);
            if (it == instance().solvers_.end()) {
                return utopia::make_unique<ConjugateGradient<Matrix, Vector>>();
            } else {
                return LinearSolverPtr(it->second->make());
            }
        }

        inline static LinearSolverPtr default_linear_solver() {
            return utopia::make_unique<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
        }

    private:
        inline static const LinearSolverFactory &instance() {
            static LinearSolverFactory instance_;
            return instance_;
        }

        LinearSolverFactory() { init(); }

        void init() {
#ifdef UTOPIA_ENABLE_LAPACK
            solvers_[Solver::direct()] = std::make_shared<LSFactoryMethod<LUDecomposition<Matrix, Vector>>>();
            solvers_[Solver::automatic()] = std::make_shared<LSFactoryMethod<LUDecomposition<Matrix, Vector>>>();
#endif  // UTOPIA_ENABLE_LAPACK
        }
    };
}  // namespace utopia

#endif  // UTOPIA_BLAS_LINEAR_SOLVER_FACTORY_HPP

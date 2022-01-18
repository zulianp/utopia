#ifndef UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP

#include "utopia_FactoryMethod.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_MatrixFreeLinearSolverFactory.hpp"
#include "utopia_Traits.hpp"

#include "utopia_petsc_Types.hpp"

#include <map>
#include <memory>
#include <string>

namespace utopia {

    template <>
    class MatrixFreeLinearSolverFactory<PetscVector, PETSC> : public BasicMatrixFreeLinearSolverFactory<PetscVector> {
    public:
        using Super = utopia::BasicMatrixFreeLinearSolverFactory<PetscVector>;
        MatrixFreeLinearSolverFactory() = default;

        static MatrixFreeLinearSolverFactory<PetscVector> &instance() {
            static MatrixFreeLinearSolverFactory<PetscVector> instance_;
            if (instance_.empty()) {
                instance_.register_solvers();
            }

            return instance_;
        }

        static MatrixFreeLinearSolverPtr new_linear_solver(const std::string &tag) {
            auto it = instance().solvers_.find(tag);
            if (it == instance().solvers_.end()) {
                return default_linear_solver();
            } else {
                return it->second->make();
            }
        }

        // static MatrixFreeLinearSolverPtr default_linear_solver() {
        //     return utopia::make_unique<ConjugateGradient<typename Traits<Vector>::Matrix, Vector, HOMEMADE>>();
        // }

        void register_solvers() override;
    };

    template <>
    class LinearSolverFactory<PetscMatrix, PetscVector, PETSC> {
    public:
        typedef utopia::LinearSolver<PetscMatrix, PetscVector> LinearSolverT;
        using LinearSolverPtr = std::unique_ptr<LinearSolverT>;
        using FactoryMethodT = utopia::IFactoryMethod<LinearSolverT>;

        template <class Alg>
        using LSFactoryMethod = FactoryMethod<LinearSolverT, Alg>;
        std::map<std::string, std::shared_ptr<FactoryMethodT>> solvers_;

        static LinearSolverPtr new_linear_solver(const std::string &tag);

        static LinearSolverPtr default_linear_solver();

        static LinearSolverFactory &instance();

    private:
        LinearSolverFactory();
        void init();
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP

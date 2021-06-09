#ifndef UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP

#include "utopia_FactoryMethod.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_Traits.hpp"

#include "utopia_petsc_Types.hpp"

#include <map>
#include <memory>
#include <string>

namespace utopia {

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

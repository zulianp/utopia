#ifndef UTOPIA_MATRIX_FREE_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_MATRIX_FREE_LINEAR_SOLVER_FACTORY_HPP

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_FactoryMethod.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_SolverType.hpp"

#include "utopia_BiCGStab.hpp"
#include "utopia_ConjugateGradient.hpp"

#include <memory>

namespace utopia {

    template <typename Vector>
    class BasicMatrixFreeLinearSolverFactory {
    public:
        using MatrixFreeLinearSolver = utopia::MatrixFreeLinearSolver<Vector>;
        using MatrixFreeLinearSolverPtr = std::unique_ptr<MatrixFreeLinearSolver>;

        using FactoryMethodT = utopia::IFactoryMethod<MatrixFreeLinearSolver>;

        template <class Alg>
        using LSFactoryMethod = FactoryMethod<MatrixFreeLinearSolver, Alg>;
        std::map<std::string, std::shared_ptr<FactoryMethodT>> solvers_;

        BasicMatrixFreeLinearSolverFactory() {}

        virtual void register_solvers() {
            register_solver<ConjugateGradient<typename Traits<Vector>::Matrix, Vector, HOMEMADE>>("cg");
            register_solver<BiCGStab<typename Traits<Vector>::Matrix, Vector, HOMEMADE>>("bcgs");
        }

        bool empty() const { return solvers_.empty(); }

        static MatrixFreeLinearSolverPtr new_linear_solver(const std::string &tag) {
            auto it = instance().solvers_.find(tag);
            if (it == instance().solvers_.end()) {
                return default_linear_solver();
            } else {
                return it->second->make();
            }
        }

        template <class Solver>
        void register_solver(const std::string &tag) {
            solvers_[tag] = utopia::make_unique<LSFactoryMethod<Solver>>();
        }

        static MatrixFreeLinearSolverPtr default_linear_solver() {
            return utopia::make_unique<ConjugateGradient<typename Traits<Vector>::Matrix, Vector, HOMEMADE>>();
        }

        static BasicMatrixFreeLinearSolverFactory<Vector> &instance() {
            static BasicMatrixFreeLinearSolverFactory<Vector> instance_;
            if (instance_.empty()) {
                instance_.register_solvers();
            }

            return instance_;
        }
    };

    /*
     * @brief      Front-end to create linear solver objects.
     *
     * @tparam     Vector
     * @tparam     Backend
     */
    template <typename Vector, int Backend = Traits<Vector>::Backend>
    class MatrixFreeLinearSolverFactory {};

    /**
     * @brief      Returns linear solver based on tag.
     *
     * @param[in]  tag     The tag - takes info on choice of LS.
     *
     * @tparam     Vector
     *
     * @return
     */
    template <class Vector>
    typename MatrixFreeLinearSolverFactory<Vector>::LinearSolverPtr matrix_free_linear_solver(
        const SolverType &tag = Solver::automatic()) {
        return MatrixFreeLinearSolverFactory<Vector>::new_linear_solver(tag);
    }

    /** \addtogroup Linear
     * @brief Solve functions for linear systems.
     * @ingroup solving
     *  @{
     */

    /**
     * @brief      Solves the linear system A * x = rhs. If no params are provided,
     * solver is choosen and set-up by default.
     *
     * @param[in]  A       The operator.
     * @param[in]  rhs     The right hand side.
     * @param      x       The initial guess/ solution.
     * @param[in]  params  The parameters.
     *
     * @tparam     Vector
     *
     * @return
     */
    template <class Vector>
    bool solve(const Operator<Vector> &A, const Vector rhs, Vector &x, const SolverType &tag = Solver::automatic()) {
        auto solver = matrix_free_linear_solver<Vector>(tag);
        solver->solve(A, rhs, x);

        return true;
    }

    /** @}*/
}  // namespace utopia

#endif  // UTOPIA_MATRIX_FREE_LINEAR_SOLVER_FACTORY_HPP
#ifndef UTOPIA_BELOS_SOLVERS_HPP
#define UTOPIA_BELOS_SOLVERS_HPP

#include "Belos_config.h"

#ifdef HAVE_BELOS_TPETRA

#include "utopia_PreconditionedSolver.hpp"
#include "utopia_trilinos_LinearSolverFactory.hpp"
#include "utopia_Smoother.hpp"

namespace utopia {
    /**@ingroup     Linear
     * @brief       Class provides interface to Trilinos Belos solvers \n
     *              For setting up basic parameters, one can use classic Belos
     * runtime options
     */
    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class BelosSolver {};
    
    template <typename Matrix, typename Vector>
    class BelosSolver<Matrix, Vector, TRILINOS> final
    : public PreconditionedSolver<Matrix, Vector>, public Smoother<Matrix, Vector> {        
    public:
        typedef UTOPIA_SCALAR(Vector) Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Preconditioner<Vector> Preconditioner;
        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;
            
        BelosSolver();
        BelosSolver(const BelosSolver &other);
        BelosSolver(Parameters params);
        ~BelosSolver();
        
        void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec) override;
        void update(const std::shared_ptr<const Matrix> &op) override;
        bool apply(const Vector &rhs, Vector &lhs) override;
        
        void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override;
        void set_preconditioner(const Matrix &precond);
        
        int get_num_iter() const;
        double achieved_tol() const;
        
        /**
         * @brief      Sets the parameters.
         *
         * @param[in]  params  The parameters
         */
        void set_parameters(const Parameters params) override;
        BelosSolver * clone() const override;
        bool smooth(const Vector &rhs, Vector &x) override;

        private:
            
            class Impl;
            std::unique_ptr<Impl> impl_;

            bool set_problem();
            bool set_problem(Matrix &A);
            void set_preconditioner();
    };

}  // namespace utopia

#endif //HAVE_BELOS_TPETRA
#endif //UTOPIA_BELOS_SOLVERS_HPP

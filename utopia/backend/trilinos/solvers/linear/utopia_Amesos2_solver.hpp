#ifndef UTOPIA_AMESOS2_SOLVERS_HPP
#define UTOPIA_AMESOS2_SOLVERS_HPP

#include "Amesos2_config.h"

#ifdef HAVE_AMESOS2_KOKKOS

#include "utopia_PreconditionedSolver.hpp"
#include "utopia_trilinos_LinearSolverFactory.hpp"
#include "utopia_Smoother.hpp"

namespace utopia {
    /**@ingroup     Linear
     * @brief       Class provides interface to Trilinos Amesos2 solvers \n
     *              For setting up basic parameters, one can use classic Amesos2
     * runtime options
     */
    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class Amesos2Solver {};
    
    template <typename Matrix, typename Vector>
    class Amesos2Solver<Matrix, Vector, TRILINOS> final
    : public PreconditionedSolver<Matrix, Vector>, public Smoother<Matrix, Vector> {        
    public:
        typedef UTOPIA_SCALAR(Vector) Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Preconditioner<Vector> Preconditioner;
        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;
            
        Amesos2Solver();
        Amesos2Solver(const Amesos2Solver &other);
        Amesos2Solver(Parameters params);
        ~Amesos2Solver();
        
        void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec) override;
        void update(const std::shared_ptr<const Matrix> &op) override;
        bool apply(const Vector &rhs, Vector &lhs) override;
        
        void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override;
        void set_preconditioner(const Matrix &precond);
        // int num_factorization() const;
        // int sym_factorization() const;
        inline int get_nnzLU () const;
        inline int get_num_preorder () const;
        inline int get_num_sym_fact () const;
        inline int get_num_numeric_fact () const;
        inline int get_num_solve () const;
        bool preordering_done () const;
        bool sym_factorization_done () const;
        bool num_factorization_done () const;

        /**
         * @brief      Checks the parameters.
         *
         * @param[in]  params  The parameters
         */
        void set_parameters(const Parameters params); // override;


         /**
          * @brief      Sets the parameters.
         */
        void check_parameters(); //override;

        Amesos2Solver * clone() const override;
        bool smooth(const Vector &rhs, Vector &x) override;

        private:
            
            class Impl;
            std::unique_ptr<Impl> impl_;

            bool set_problem();
            bool set_problem(Matrix &A);
            void set_preconditioner();
    };

}  // namespace utopia

#endif //HAVE_AMESOS2_KOKKOS
#endif //UTOPIA_AMESOS2_SOLVERS_HPP

#ifndef UTOPIA_BELOS_SOLVERS_HPP
#define UTOPIA_BELOS_SOLVERS_HPP

#include <utopia_Config.hpp>

#ifdef WITH_TRILINOS_BELOS

#include "Belos_config.h"

#include "utopia_PreconditionedSolver.hpp"
#include "utopia_trilinos_LinearSolverFactory.hpp"
#include "utopia_Smoother.hpp"

namespace utopia {
    /**@ingroup     Linear
     * @brief       Class provides interface to Trilinos Belos solvers \n
     *              For setting up basic parameters, one can use classic Belos
     *              runtime options
     */
    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class BelosSolver {};

    template <typename Matrix, typename Vector>
    class BelosSolver<Matrix, Vector, TRILINOS> final : public PreconditionedSolver<Matrix, Vector>
    {
    public:
        typedef UTOPIA_SCALAR(Vector) Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        
        typedef utopia::Preconditioner<Vector> Preconditioner;
        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;

        BelosSolver();
        BelosSolver(const BelosSolver &other);
        ~BelosSolver() override;

        void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec) override;
        void update(const std::shared_ptr<const Matrix> &op) override;
        bool apply(const Vector &rhs, Vector &lhs) override;

        void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override;
        void set_preconditioner(const Matrix &precond);

        int get_num_iter() const;
        double achieved_tol() const;

        /**
         * @brief      Reads the xml file based on different layout than read
         *
         * @param[path] path  location of the xml file
         */
        void read_xml(const std::string &path);

        /**
         * @brief      Reads the UTOPIA input object.
         *
         * @param[in]  in  The parameters
         */
        void read(Input &in) override;/* {
          Smoother<Matrix, Vector>::read(in); 
          PreconditionedSolver::read(in);
          //TODO
          atol_(1e-9);
          rtol_(1e-9),;
          stol_(1e-11);
          max_it_(300);
          verbose_(false);

          m_utopia_warning_once("not implemented");
        }*/

        /**
         * @brief      Prints the parameters used.
         *
         * @param[os]  os The std::ostream parameters
         */
        void print_usage(std::ostream &os = std::cout) const override;
   /*     {
          Smoother<Matrix, Vector>::print_usage(os); 
          PreconditionedSolver::print_usage(os);
          //TODO
          m_utopia_warning_once("not implemented");

        }*/

        /**
         * @brief      Clone the object.
         *
         */
        BelosSolver * clone() const override;
        

        private:

            class Impl;
            std::unique_ptr<Impl> impl_;

            bool set_problem();
            bool set_problem(Matrix &A);
            void set_preconditioner();
    };

}  // namespace utopia

#else  // WITH_TRILINOS_BELOS
  #error "Trilinos was not configured with Belos, hence you cannot use the utopia::BelosSolver."
#endif //WITH_TRILINOS_BELOS
#endif //UTOPIA_BELOS_SOLVERS_HPP

#ifndef UTOPIA_TRILINOS_SOLVERS_HPP
#define UTOPIA_TRILINOS_SOLVERS_HPP

//TODO here we are goint to put the relevant hpp for trilinos solvers
#include "utopia_trilinos_LinearSolverFactory.hpp"


namespace utopia
{


/**@ingroup     Linear
 * @brief       Class provides interface to Trilinos Belos solvers \n
 *              For setting up basic parameters, one can use classic Belos runtime options
 */
template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
class BelosSolver {};


template<typename Matrix, typename Vector>
class BelosSolver<Matrix, Vector, Trilinos> : virtual public PreconditionedSolver<Matrix, Vector>
    {

    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Preconditioner<Vector> Preconditioner;
        // typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;


        virtual ~BelosSolver()
            {
            }


        /**
         * @brief      Sets the parameters.
         *
         * @param[in]  params  The parameters
         */
        virtual void set_parameters(const Parameters params) override
            {
            PreconditionedSolver::set_parameters(params);
            }

    }

#endif //UTOPIA_TRILINOS_SOLVERS_HPP

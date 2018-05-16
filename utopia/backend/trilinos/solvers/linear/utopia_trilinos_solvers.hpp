#ifndef UTOPIA_TRILINOS_SOLVERS_HPP
#define UTOPIA_TRILINOS_SOLVERS_HPP

//TODO here we are goint to put the relevant hpp for trilinos solvers
#include "utopia_trilinos_LinearSolverFactory.hpp"
#include "utopia_PreconditionedSolver.hpp"

//#include "BelosConfigDefs.hpp"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
//#include "BelosBlockCGSolMgr.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_DefaultNode.hpp>


namespace utopia
{

         typedef double                           ST;
         typedef Teuchos::ScalarTraits<ST>        SCT;
         typedef SCT::magnitudeType               MT;
         typedef Tpetra::Operator<ST,int>         OP;
         typedef Tpetra::MultiVector<ST,int>      MV;
         typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
         typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::LinearProblem<SC, MV, OP> problem_type;
typedef Belos::SolverManager<SC, MV, OP> solver_type;

/**@ingroup     Linear
 * @brief       Class provides interface to Trilinos Belos solvers \n
 *              For setting up basic parameters, one can use classic Belos runtime options
 */
    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class BelosSolver {};


    template<typename Matrix, typename Vector>
    class BelosSolver<Matrix, Vector, TRILINOS> : virtual public PreconditionedSolver<Matrix, Vector>
    {

//        Belos::LinearProblem<ST,MV,OP> problem;
Teuchos::RCP<Belos::LinearProblem<SC, MV, OP> > linearProblem;

//
Teuchos::RCP<Teuchos::ParameterList> ParamList;
Teuchos::RCP<solver_type> belosSolver;
Belos::SolverFactory<SC, MV, OP> belosFactory;
        bool success = false;
        bool verbose = false;


    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Preconditioner<Vector> Preconditioner;
        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;

//////
BelosSolver(A,LHS,RHS, param_file_name)
{
  linearProblem = 
 Teuchos::rcp (new problem_type (A, LHS, RHS ));
ParamList = Teuchos::getParametersFromXmlFile(param_file_name);
}

setProblem(it_sol_type,ParamList)
{
linearProblem->setProblem ();
 belosSolver = belosFactory.create (it_sol_type, ParamList);
//linearProblem->setRightPrec (M_ifpack);
belosSolver->setProblem (linearProblem) ;
}

getNumIter(int numIter)
{
numIter = belosSolver->getNumIters();
}
/////
        virtual ~BelosSolver()
        {
        }


        bool apply(const Vector &rhs, Vector &sol) override
	{
            linearProblem->apply(*rhs, *sol);
            belosSolver->solve ();
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

    };
}

#endif //UTOPIA_TRILINOS_SOLVERS_HPP

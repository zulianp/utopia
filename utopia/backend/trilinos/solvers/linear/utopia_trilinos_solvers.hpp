#ifndef UTOPIA_TRILINOS_SOLVERS_HPP
#define UTOPIA_TRILINOS_SOLVERS_HPP

// TODO here we are goint to put the relevant hpp for trilinos solvers
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_trilinos_LinearSolverFactory.hpp"

//#include "BelosConfigDefs.hpp"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
//#include "BelosBlockCGSolMgr.hpp"

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>

namespace utopia {

typedef double ST;
typedef Teuchos::ScalarTraits<ST> SCT;
typedef SCT::magnitudeType MT;

typedef Tpetra::Operator<ST, int> OP;
typedef Tpetra::MultiVector<ST, int> MV;
typedef Tpetra::Operator<>::scalar_type SC;
typedef Tpetra::Operator<SC>::local_ordinal_type LO;
typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;

typedef Belos::OperatorTraits<ST, MV, OP> OPT;
typedef Belos::MultiVecTraits<ST, MV> MVT;
typedef Belos::LinearProblem<SC, MV, OP> problem_type;
typedef Belos::SolverManager<SC, MV, OP> solver_type;

typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;
typedef Kokkos::Compat::KokkosThreadsWrapperNode thread_node;

typedef serial_node NT;

typedef Tpetra::Vector<SC, LO, GO, NT> vec_type;
typedef Tpetra::CrsMatrix<SC, LO, GO, NT> matrix_type;

/**@ingroup     Linear
 * @brief       Class provides interface to Trilinos Belos solvers \n
 *              For setting up basic parameters, one can use classic Belos
 * runtime options
 */
template <typename Matrix, typename Vector,
          int Backend = Traits<Matrix>::Backend>
class BelosSolver {};

template <typename Matrix, typename Vector>
class BelosSolver<Matrix, Vector, TRILINOS>
    : virtual public PreconditionedSolver<Matrix, Vector> {
  //        Belos::LinearProblem<ST,MV,OP> problem;
  Teuchos::RCP<Belos::LinearProblem<SC, MV, OP> > linearProblem;
  //
  Teuchos::RCP<Teuchos::ParameterList> ParamList;
  Teuchos::RCP<solver_type> belosSolver;
  Belos::SolverFactory<SC, MV, OP> belosFactory;
  //
  bool success = false;
  bool verbose = false;

 public:
  typedef UTOPIA_SCALAR(Vector) Scalar;
  typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
  typedef utopia::Preconditioner<Vector> Preconditioner;
  typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
  typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
  typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;

  //////
  BelosSolver(Teuchos::RCP<matrix_type> A, Teuchos::RCP<vec_type> LHS,
              Teuchos::RCP<vec_type> RHS, std::string param_file_name) {
    linearProblem = Teuchos::rcp(new problem_type(A, LHS, RHS));
    ParamList = Teuchos::getParametersFromXmlFile(param_file_name);
  }

  bool setProblem(std::string it_sol_type,
                  Teuchos::RCP<Teuchos::ParameterList> ParamList) {
    linearProblem->setProblem();
    belosSolver = belosFactory.create(it_sol_type, ParamList);
    // linearProblem->setRightPrec (M_ifpack);
    belosSolver->setProblem(linearProblem);
    return true;
  }

  /////
  virtual ~BelosSolver() {}

  bool apply(const Vector &rhs, Vector &sol) override {
    linearProblem->apply(*rhs, *sol);
    belosSolver->solve();
    return true;
  }

  int getNumIter() { return belosSolver->getNumIters(); }

  /**
   * @brief      Sets the parameters.
   *
   * @param[in]  params  The parameters
   */
  virtual void set_parameters(const Parameters params,
                              std::string it_sol_type) override {
    PreconditionedSolver::set_parameters(params);
    setProblem(it_sol_type, ParamList);
  }
};
}  // namespace utopia

#endif  // UTOPIA_TRILINOS_SOLVERS_HPP

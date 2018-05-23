#ifndef UTOPIA_BELOS_SOLVERS_HPP
#define UTOPIA_BELOS_SOLVERS_HPP

#include "utopia_PreconditionedSolver.hpp"
#include "utopia_trilinos_LinearSolverFactory.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp> //TODO remove from here


#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>

//TODO remove from here
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>

/*#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
*/
#include <Ifpack2_Factory.hpp>


namespace utopia {

typedef double ST;
//typedef Teuchos::ScalarTraits<ST> SCT;
//typedef SCT::magnitudeType MT;

typedef Tpetra::Operator<>::scalar_type SC;
typedef Tpetra::Operator<SC>::local_ordinal_type LO;
typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;
//typedef Tpetra::Map<LO, GO, NT> map_type;

typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;
typedef Kokkos::Compat::KokkosThreadsWrapperNode thread_node;

typedef serial_node NT;

typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
typedef Tpetra::Operator<SC, LO, GO, NT> OP;


typedef Belos::LinearProblem<SC, MV, OP> problem_type;
typedef Belos::SolverManager<SC, MV, OP> solver_type;

typedef Ifpack2::Preconditioner<SC, LO, GO, NT> ifpack_prec_type;

//typedef MueLu::TpetraOperator<SC, LO, GO, NT> muelu_prec_type;

/*
typedef Tpetra::Vector<SC, LO, GO, NT> vec_type;*/
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
    : public PreconditionedSolver<Matrix, Vector>, public Smoother<Matrix, Vector> {

 private:
  Teuchos::RCP<problem_type> linearProblem;
  Teuchos::RCP<Teuchos::ParameterList> ParamList;
//  auto& utopiaPL;// ParamList->sublist("UTOPIA", true);
  Teuchos::RCP<solver_type> belosSolver;
  Belos::SolverFactory<SC, MV, OP> belosFactory;

  //preconditioner
  Teuchos::RCP<ifpack_prec_type> M_ifpack;
//  Teuchos::RCP<muelu_prec_type> M_muelu;

 public:
  typedef UTOPIA_SCALAR(Vector) Scalar;
  typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
  typedef utopia::Preconditioner<Vector> Preconditioner;
  typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
  typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
  typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;

  //////

  BelosSolver( Parameters params)
  {
  set_parameters(params);
  }

/*  BelosSolver(Matrix A, Vector LHS,     ///Teuchos::RCP<matrix_type> A, Teuchos::RCP<vec_type> LHS,
              Vector RHS) {
    linearProblem = Teuchos::rcp(new problem_type(A, LHS, RHS));
  }*/

  bool setProblem() {
    linearProblem->setProblem();
    belosSolver = belosFactory.create( ParamList->sublist("UTOPIA", true).get("Solver Type", "CG"), ParamList); //to change it to have the specialization
    belosSolver->setProblem(linearProblem);
    if (this->verbose()) belosSolver->getCurrentParameters()->print();

    return true;
  }

  /////
  virtual ~BelosSolver() {}

  bool apply(const Vector &RHS, Vector &LHS) override {
 if (linearProblem.is_null())
   {linearProblem = Teuchos::rcp(new problem_type(this->get_operator(), LHS, RHS)); }

   setProblem();
   belosSolver->solve();

    return true;
  }
   bool solve(Matrix &A, const Vector &RHS, Vector &LHS) //override 
{
  if (linearProblem.is_null())
     {linearProblem = Teuchos::rcp(new problem_type(A, LHS, RHS));}
     set_preconditioner(A);
     setProblem();
     belosSolver->solve();
     return true;
   } 

void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override
{
bool direct_solver = ParamList->sublist("UTOPIA", true).get<bool>("Direct Preconditioner", false);
std::string dir_prec_type = ParamList->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");

      if ( direct_solver ) {
        M_ifpack = Ifpack2::Factory::create<matrix_type>(dir_prec_type, *precond);
        assert(!M_ifpack.is_null());
        M_ifpack->setParameters(ParamList->sublist(dir_prec_type, false));
        M_ifpack->initialize();
        M_ifpack->compute();
        linearProblem->setLeftPrec(M_ifpack);
      } else {
        // Multigrid Hierarchy
        //M_muelu = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<OP>)*precond,
        //                                            ParamList->sublist("MueLu", false));
        //assert(!M_muelu.is_null());
        //linearProblem->setRightPrec(M_muelu);
      }
}

void set_preconditioner(const Matrix &precond) //override
{
bool direct_solver = ParamList->sublist("UTOPIA", true).get<bool>("Direct Preconditioner", false);
std::string dir_prec_type = ParamList->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");
      
      if ( direct_solver ) {
        M_ifpack = Ifpack2::Factory::create<matrix_type>(dir_prec_type, *precond);
        assert(!M_ifpack.is_null());
        M_ifpack->setParameters(ParamList->sublist(dir_prec_type, false));
        M_ifpack->initialize();
        M_ifpack->compute();
        linearProblem->setLeftPrec(M_ifpack);
     } else {
       // Multigrid Hierarchy
       //M_muelu = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<OP>)*precond,
       //                                            ParamList->sublist("MueLu", false));
       //assert(!M_muelu.is_null());
       //linearProblem->setRightPrec(M_muelu);
}}



  int getNumIter() { return belosSolver->getNumIters(); }
  double achievedTol() { return belosSolver->achievedTol(); }

  /**
   * @brief      Sets the parameters.
   *
   * @param[in]  params  The parameters
   */
   void set_parameters(const Parameters params){
    ParamList = Teuchos::getParametersFromXmlFile(params.param_file_name);
  }

                     virtual BelosSolver * clone() const override
                     {
                        return new BelosSolver(*this);
                     }

                      void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec) override
                     {
                       PreconditionedSolver::update(op, prec); 
                     }

                   void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Vector> &prec) override
{
PreconditionedSolver::update(op, prec);
}

                     bool smooth(const Vector &rhs, Vector &x) override
                     {
                        return false;
                     }


};
}  // namespace utopia

#endif  // UTOPIA_BELOS_SOLVERS_HPP

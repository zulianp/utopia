#ifndef UTOPIA_BELOS_SOLVERS_HPP
#define UTOPIA_BELOS_SOLVERS_HPP

#include "Belos_config.h"

#ifdef HAVE_BELOS_TPETRA

#include "utopia_PreconditionedSolver.hpp"
#include "utopia_trilinos_LinearSolverFactory.hpp"
#include "utopia_Smoother.hpp"

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


//FIXME find right macros (these packages are optional in trilinos, they should be optional also in utopia)
#define HAVE_BELOS_MUELU
#define HAVE_BELOS_IFPACK2


#ifdef HAVE_BELOS_MUELU
#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
#else
#warning "HAVE_BELOS_MUELU not defined"
#endif //HAVE_BELOS_MUELU


#ifdef HAVE_BELOS_IFPACK2
#include <Ifpack2_Factory.hpp>
#else
#warning "HAVE_BELOS_IFPACK2 not defined"
#endif //HAVE_BELOS_IFPACK


namespace utopia {
    
    typedef double ST;
    //typedef Teuchos::ScalarTraits<ST> SCT;
    //typedef SCT::magnitudeType MT;
    
    typedef Tpetra::Operator<>::scalar_type SC;
    typedef Tpetra::Operator<SC>::local_ordinal_type LO;
    typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;
    //typedef Tpetra::Map<LO, GO, NT> map_type;
    
    typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;
    
#ifdef  KOKKOS_CUDA
    typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
    typedef cuda_node NT;
#elif defined   KOKKOS_OPENMP
    typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
    typedef Kokkos::Compat::KokkosThreadsWrapperNode thread_node;
    typedef openmp_node NT;
#else
    typedef serial_node NT;
#endif
    
    typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
    typedef Tpetra::Operator<SC, LO, GO, NT> OP;
    
    
    typedef Belos::LinearProblem<SC, MV, OP> problem_type;
    typedef Belos::SolverManager<SC, MV, OP> solver_type;
    
#ifdef HAVE_BELOS_IFPACK2
    typedef Ifpack2::Preconditioner<SC, LO, GO, NT> ifpack_prec_type;
#endif //HAVE_BELOS_IFPACK
    
#ifdef HAVE_BELOS_MUELU
    typedef MueLu::TpetraOperator<SC, LO, GO, NT> muelu_prec_type;
#endif
    
    
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
    class BelosSolver<Matrix, Vector, TRILINOS> final
    : public PreconditionedSolver<Matrix, Vector>, public Smoother<Matrix, Vector> {
        
    private:
        Teuchos::RCP<problem_type> linearProblem;
        Teuchos::RCP<Teuchos::ParameterList> ParamList;
        //  auto& utopiaPL;// ParamList->sublist("UTOPIA", true);
        Teuchos::RCP<solver_type> belosSolver;
        Belos::SolverFactory<SC, MV, OP> belosFactory;
        
        //preconditioner
#ifdef HAVE_BELOS_IFPACK2
        Teuchos::RCP<ifpack_prec_type> M_ifpack;
#endif //HAVE_BELOS_IFPACK2
        
#ifdef HAVE_BELOS_MUELU
        Teuchos::RCP<muelu_prec_type> M_muelu;
#endif //HAVE_BELOS_MUELU
        
    public:
        typedef UTOPIA_SCALAR(Vector) Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Preconditioner<Vector> Preconditioner;
        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;
        
        //////
        BelosSolver()
        {

        }

        BelosSolver(Parameters params)
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

        bool setProblem(Matrix &A) {
             linearProblem->setProblem();
            belosSolver = belosFactory.create( ParamList->sublist("UTOPIA", true).get("Solver Type", "CG"), ParamList); //to change it to have the specialization
            set_preconditioner(A);
            belosSolver->setProblem(linearProblem);
             if (this->verbose()) belosSolver->getCurrentParameters()->print();

             return true;
         }

        /////
        virtual ~BelosSolver() {}
        
        
        void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec) override
        {
            PreconditionedSolver::update(op, prec);
            //TODO
        }
        
        void update(const std::shared_ptr<const Matrix> &op) override
        {
            PreconditionedSolver::update(op);
            //TODO
        }
        
        
        bool apply(const Vector &rhs, Vector &lhs) override {
            linearProblem = Teuchos::rcp(
             new problem_type(
               this->get_operator()->implementation().implementation_ptr(),
               lhs.implementation().implementation_ptr(),
               rhs.implementation().implementation_ptr()
               )
             );
            setProblem();
            belosSolver->solve();
            return true;
        }
        
bool solve(Matrix &A,const Vector &rhs, Vector &lhs) { //override {
        linearProblem = Teuchos::rcp( new problem_type(A.implementation().implementation_ptr(),
        lhs.implementation().implementation_ptr(),
        rhs.implementation().implementation_ptr() ) );
            setProblem(A);
            belosSolver->solve();
            return true;
        }
        
        
        void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override
        {
            bool direct_solver = ParamList->sublist("UTOPIA", true).get<bool>("Direct Preconditioner", false);
            std::string dir_prec_type = ParamList->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");
            //TODO    
            auto delegate_ptr = std::dynamic_pointer_cast<DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());
            /*if(delegate_ptr) {
             if(ksp_->has_shell_pc()) {
             m_utopia_warning_once("set_preconditioner sets jacobi if a delegate precond has been set and type is matshell");
             ksp_->pc_type("jacobi");
             }                                                                                                                                                                                                                  
             } else if(this->get_preconditioner()) {
             auto shell_ptr = this->get_preconditioner().get();
             ksp_->attach_shell_preconditioner(UtopiaPCApplyShell,
             shell_ptr,
             nullptr,
             nullptr
             );
             }
             }*/
            
            if ( direct_solver ) 
            {
                //                 M_ifpack = Ifpack2::Factory::create<matrix_type>(dir_prec_type, precond->implementation().implementation_ptr()); //TODO
                
#ifdef HAVE_BELOS_IFPACK2
                assert(!M_ifpack.is_null());
                M_ifpack->setParameters(ParamList->sublist(dir_prec_type, false));
                M_ifpack->initialize();
                M_ifpack->compute();
                linearProblem->setLeftPrec(M_ifpack);
#endif //HAVE_BELOS_IFPACK2
                
            } else {
#ifdef HAVE_BELOS_MUELU
                // Multigrid Hierarchy
                //                M_muelu = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<OP>)precond->implementation().implementation_ptr(),    //TODO
                //                                                            ParamList->sublist("MueLu", false));
                assert(!M_muelu.is_null());
                linearProblem->setRightPrec(M_muelu);
#else
                assert(false);
#endif //HAVE_BELOS_MUELU
            }
        }
        
        void set_preconditioner(const Matrix &precond) //override
        {
            bool direct_solver = ParamList->sublist("UTOPIA", true).get<bool>("Direct Preconditioner", false);
            std::string dir_prec_type = ParamList->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");
            if ( direct_solver ) {
#ifdef HAVE_BELOS_IFPACK2
                M_ifpack = Ifpack2::Factory::create<matrix_type>(dir_prec_type, precond.implementation().implementation_ptr());
                assert(!M_ifpack.is_null());
                M_ifpack->setParameters(ParamList->sublist(dir_prec_type, false));
                M_ifpack->initialize();
                M_ifpack->compute();
                linearProblem->setLeftPrec(M_ifpack);
#endif //HAVE_BELOS_IFPACK2
            } else {
                
#ifdef HAVE_BELOS_MUELU
                // Multigrid Hierarchy
                M_muelu = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<OP>) precond.implementation().implementation_ptr(),
                                                            ParamList->sublist("MueLu", false));
                assert(!M_muelu.is_null());
                linearProblem->setRightPrec(M_muelu);
#else
                assert(false);
#endif //HAVE_BELOS_MUELU
            }
        }
        
        
        int getNumIter() { return belosSolver->getNumIters(); }
        double achievedTol() { return belosSolver->achievedTol(); }
        
        /**
         * @brief      Sets the parameters.
         *
         * @param[in]  params  The parameters
         */
        void set_parameters(const Parameters params) override {
            if(!params.param_file_name().empty()) {
                try {
                    ParamList = Teuchos::getParametersFromXmlFile(params.param_file_name());
                } catch(const std::exception &ex) {
                    std::cerr << ex.what() << std::endl;
                    assert(false);
                    abort();
                }
            } else {
                //use default paramlist
            }
        }
        
        virtual BelosSolver * clone() const override
        {
            return new BelosSolver(*this);
        }
        
        
        bool smooth(const Vector &rhs, Vector &x) override
        {
            return false;
        }
    };
}  // namespace utopia

#endif //HAVE_BELOS_TPETRA
#endif //UTOPIA_BELOS_SOLVERS_HPP

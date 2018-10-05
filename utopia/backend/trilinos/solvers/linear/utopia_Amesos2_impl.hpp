#ifndef UTOPIA_AMESOS2_IMPL_HPP
#define UTOPIA_AMESOS2_IMPL_HPP

#include "utopia_Amesos2_solver.hpp"

#include "utopia_make_unique.hpp"

#include <Amesos2_Factory.hpp>
#include <Amesos2_Solver.hpp>

//TODO remove from here
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>


// #ifdef HAVE_AMESOS2_TPETRA

//FIXME find right macros (these packages are optional in trilinos, they should be optional also in utopia)
#define HAVE_AMESOS2_MUELU
#define HAVE_AMESOS2_IFPACK2


// #ifdef HAVE_AMESOS2_MUELU
// #include <MueLu.hpp>
// #include <MueLu_CreateTpetraPreconditioner.hpp>
// #include <MueLu_TpetraOperator.hpp>
// #else
// #warning "HAVE_AMESOS2_MUELU not defined"
// #endif //HAVE_AMESOS2_MUELU


// #ifdef HAVE_AMESOS2_IFPACK2
// #include <Ifpack2_Factory.hpp>
// #else
// #warning "HAVE_AMESOS2_IFPACK2 not defined"
// #endif //HAVE_AMESOS2_IFPACK


namespace utopia {
    
    template <typename Matrix, typename Vector>
    class Amesos2Solver<Matrix, Vector, TRILINOS>::Impl {
    public:
        typedef double ST;
        
        typedef Tpetra::Operator<>::scalar_type SC;
        typedef Tpetra::Operator<SC>::local_ordinal_type LO;
        typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;
        
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
        
        typedef Amesos2::LinearProblem<SC, MV, OP> problem_type;
        typedef Amesos2::SolverManager<SC, MV, OP> solver_type;
        
// #ifdef HAVE_AMESOS2_IFPACK2
        typedef Ifpack2::Preconditioner<SC, LO, GO, NT> ifpack_prec_type;
// #endif //HAVE_AMESOS2_IFPACK
        
// #ifdef HAVE_AMESOS2_MUELU
        typedef MueLu::TpetraOperator<SC, LO, GO, NT> muelu_prec_type;
// #endif
        
        typedef Tpetra::Vector<SC, LO, GO, NT> vec_type;
        typedef Tpetra::CrsMatrix<SC, LO, GO, NT> matrix_type;
        
        Teuchos::RCP<problem_type> linear_problem;
        Teuchos::RCP<Teuchos::ParameterList> param_list;
        //  auto& utopiaPL;// impl_->param_list->sublist("UTOPIA", true);
        Teuchos::RCP<solver_type> amesos2_solver;
        Amesos2::SolverFactory<SC, MV, OP> amesos2_factory;
        
        //preconditioner
// #ifdef HAVE_AMESOS2_IFPACK2
        Teuchos::RCP<ifpack_prec_type> M_ifpack;
// #endif //HAVE_AMESOS2_IFPACK2
        
// #ifdef HAVE_AMESOS2_MUELU
        Teuchos::RCP<muelu_prec_type> M_muelu;
// #endif //HAVE_AMESOS2_MUELU
        
    };
    
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::Amesos2Solver(Parameters params)
    : impl_(make_unique<Impl>())
    {
        set_parameters(params);
    }

    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::Amesos2Solver(const Amesos2Solver &other)
    : impl_(make_unique<Impl>(*other.impl_)) {
        //FIXME
    }
    
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::~Amesos2Solver() {}
    
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::Amesos2Solver() : impl_(make_unique<Impl>()) {}
    
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op,
                                                       const std::shared_ptr<const Matrix> &prec)
    {
        PreconditionedSolver::update(op, prec);
        // set_problem(*op);
    }
    
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op)
    {
        PreconditionedSolver::update(op);
        // set_problem(*op);
    }
    
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::apply(const Vector &rhs, Vector &lhs) {
        
        impl_->linear_problem = Teuchos::rcp(
                                             new typename Impl::problem_type(
                                                                             this->get_operator()->implementation().implementation_ptr(),
                                                                             lhs.implementation().implementation_ptr(),
                                                                             rhs.implementation().implementation_ptr()
                                                                             )
                                             );

        set_problem();

        assert((impl_->amesos2_solver));
        impl_->amesos2_solver->solve();
        return true;
    }
    
    template <typename Matrix, typename Vector>
    int Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_iter() const { return impl_->amesos2_solver->getNumIters(); }
    
    
    template <typename Matrix, typename Vector>
    double Amesos2Solver<Matrix, Vector, TRILINOS>::achieved_tol() const { return impl_->amesos2_solver->achievedTol(); }
    
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
    {
        set_preconditioner(); //(A) //FIXME
    }
    
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::set_preconditioner(const Matrix &precond)
    {
        bool direct_solver = impl_->param_list->sublist("UTOPIA", true).template get<bool>("Direct Preconditioner", false);
        std::string dir_prec_type = impl_->param_list->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");
        if ( direct_solver ) {
// #ifdef HAVE_AMESOS2_IFPACK2
            impl_->M_ifpack = Ifpack2::Factory::create<typename Impl::matrix_type>(dir_prec_type, precond.implementation().implementation_ptr());
            assert(!impl_->M_ifpack.is_null());
            impl_->M_ifpack->setParameters(impl_->param_list->sublist(dir_prec_type, false));
            impl_->M_ifpack->initialize();
            impl_->M_ifpack->compute();
            impl_->linear_problem->setLeftPrec(impl_->M_ifpack);
// #endif //HAVE_AMESOS2_IFPACK2
        } else {
// #ifdef HAVE_AMESOS2_MUELU
            // Multigrid Hierarchy
            impl_->M_muelu = MueLu::CreateTpetraPreconditioner((
                                                                Teuchos::RCP<typename Impl::OP>) precond.implementation().implementation_ptr(),
                                                               impl_->param_list->sublist("MueLu", false)
                                                               );
            
            assert(!impl_->M_muelu.is_null());
            impl_->linear_problem->setRightPrec(impl_->M_muelu);
// #else
            assert(false);
// #endif //HAVE_AMESOS2_MUELU
        }
    }
    
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::set_parameters(const Parameters params)
    {
        if(!params.param_file_name().empty()) {
            try {
                impl_->param_list = Teuchos::getParametersFromXmlFile(params.param_file_name());
            } catch(const std::exception &ex) {
                std::cerr << ex.what() << std::endl;
                assert(false);
                abort();
            }
        } else {
            //use default paramlist
        }
    }
    
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS> * Amesos2Solver<Matrix, Vector, TRILINOS>::clone() const
    {
        return new Amesos2Solver(*this);
    }
    
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::smooth(const Vector &rhs, Vector &x)
    {
        return false;
    }
    
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::set_problem()
    {
        impl_->linear_problem->setProblem();
        impl_->amesos2_solver = impl_->amesos2_factory.create( impl_->param_list->sublist("UTOPIA", true).get("Solver Type", "CG"), impl_->param_list); //to change it to have the specialization
        impl_->amesos2_solver->setProblem(impl_->linear_problem);
        if (this->verbose()) { impl_->amesos2_solver->getCurrentParameters()->print(); }
        
        return true;
    }
    
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::set_problem(Matrix &A)
    {
        impl_->linear_problem->setProblem();
        impl_->amesos2_solver = impl_->amesos2_factory.create( impl_->param_list->sublist("UTOPIA", true).get("Solver Type", "CG"), impl_->param_list); //to change it to have the specialization
        set_preconditioner(); //(A);
        impl_->amesos2_solver->setProblem(impl_->linear_problem);
        if (this->verbose()) { impl_->amesos2_solver->getCurrentParameters()->print(); }
        return true;
    }
    
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::set_preconditioner()//const std::shared_ptr<Preconditioner> &precond)
    {
        bool direct_solver = impl_->param_list->sublist("UTOPIA", true).template get<bool>("Direct Preconditioner", false);
        std::string dir_prec_type = impl_->param_list->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");
        
        if ( direct_solver )
        {
            //                 M_ifpack = Ifpack2::Factory::create<matrix_type>(dir_prec_type, precond->implementation().implementation_ptr()); //TODO
            
// #ifdef HAVE_AMESOS2_IFPACK2
            assert(!impl_->M_ifpack.is_null());
            impl_->M_ifpack->setParameters(impl_->param_list->sublist(dir_prec_type, false));
            impl_->M_ifpack->initialize();
            impl_->M_ifpack->compute();
            impl_->linear_problem->setLeftPrec(impl_->M_ifpack);
// #endif //HAVE_AMESOS2_IFPACK2
            
        } else {
// #ifdef HAVE_AMESOS2_MUELU
            // Multigrid Hierarchy
            //                M_muelu = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<OP>)precond->implementation().implementation_ptr(),    //TODO
            //                                                            impl_->param_list->sublist("MueLu", false));
            assert(!impl_->M_muelu.is_null());
            impl_->linear_problem->setRightPrec(impl_->M_muelu);
// #else
            assert(false);
// #endif //HAVE_AMESOS2_MUELU
        }
    }
    
}  // namespace utopia

// #endif //HAVE_AMESOS2_TPETRA
#endif //UTOPIA_AMESOS2_IMPL_HPP

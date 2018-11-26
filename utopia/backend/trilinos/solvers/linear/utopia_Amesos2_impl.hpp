#ifndef UTOPIA_AMESOS2_IMPL_HPP
#define UTOPIA_AMESOS2_IMPL_HPP

#include "utopia_Amesos2_solver.hpp"

#include "utopia_make_unique.hpp"

#include "Amesos2.hpp"
#include "Amesos2_Meta.hpp"

//TODO remove from here
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

////
#include <Amesos2_Factory.hpp>
#include <Amesos2_Solver.hpp>
#include <Amesos2_MultiVecAdapter.hpp>

#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverManager.hpp>

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Kokkos_Macros.hpp>
#include <Kokkos_DefaultNode.hpp>

#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>

//////



#ifdef HAVE_AMESOS2_KOKKOS

//FIXME find right macros (these packages are optional in trilinos, they should be optional also in utopia)
// #define HAVE_AMESOS2_MUELU
// #define HAVE_AMESOS2_IFPACK2


#ifdef WITH_TRILINOS_MUELU
#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_Utilities.hpp>
#else
#warning "WITH_TRILINOS_MUELU not defined"
#endif //WITH_TRILINOS_MUELU


#ifdef WITH_TRILINOS_IFPACK2
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Preconditioner.hpp>
#else
#warning "WITH_TRILINOS_IFPACK2 not defined"
#endif //WITH_TRILINOS_IFPACK2


namespace utopia {
    /**
     * The Amesos2Solver class is the implementaion class
     *
     * \author Nur Aiman Fadel
     */
    template <typename Matrix, typename Vector>
    class Amesos2Solver<Matrix, Vector, TRILINOS>::Impl {
    public:

      //FIXME change all these typedefs by accessing inner definitions of matrix and vector
        //example
        using vec_impl          = typename Vector::Implementation;
        using multi_vector_type = typename vec_impl::multi_vector_type;

        //....

        typedef double ST;
        
        typedef Tpetra::Operator<>::scalar_type SC;
        typedef Tpetra::Operator<>::local_ordinal_type LO;
        typedef Tpetra::Operator<>::global_ordinal_type GO;
        
        typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;
        
#ifdef  KOKKOS_CUDA
        typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
        typedef cuda_node NT;
#elif defined  KOKKOS_ROCM //Kokkos::Compat::KokkosROCmWrapperNode doesn't exist
        typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::ROCm> rocm_node;
        typedef rocm_node NT;
#elif defined   KOKKOS_OPENMP
        typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
        typedef Kokkos::Compat::KokkosThreadsWrapperNode thread_node;
        typedef openmp_node NT;
#else
        typedef serial_node NT;
#endif
        
        typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
        typedef Tpetra::Operator<SC, LO, GO, NT> OP;
        
        // typedef Tpetra::Vector<SC, LO, GO, NT> vec_type;
        
        typedef Tpetra::CrsMatrix<SC, LO, GO, NT> matrix_type;
        
        //typedef Amesos2::LinearProblem<SC, MV, OP> problem_type;
        typedef Amesos2::Solver<matrix_type, multi_vector_type> solver_type;
        #ifdef WITH_TRILINOS_IFPACK2
        typedef Ifpack2::Preconditioner<SC, LO, GO, NT> ifpack_prec_type;
        #endif //WITH_TRILINOS_IFPACK2
        
        #ifdef WITH_TRILINOS_MUELU
        typedef MueLu::TpetraOperator<SC, LO, GO, NT> muelu_prec_type;
        #endif //WITH_TRILINOS_MUELU
        
        // Members
        //Teuchos::RCP<problem_type> linear_problem;
        Teuchos::RCP<Teuchos::ParameterList> param_list_;
        //  auto& utopiaPL;// impl_->param_list_->sublist("UTOPIA", true);
        Teuchos::RCP<solver_type> solver_;
        
        //preconditioner
        #ifdef WITH_TRILINOS_IFPACK2
        Teuchos::RCP<ifpack_prec_type> ifpack_prec_;
        #endif //WITH_TRILINOS_IFPACK2
        
        #ifdef WITH_TRILINOS_MUELU
        Teuchos::RCP<muelu_prec_type> muelu_prec_;
        #endif //WITH_TRILINOS_MUELU
        
    };
    
    /**
     * Constructor that sets the solver's parameters.
     * \param params an object Parameters
     */
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::Amesos2Solver(Parameters params)
    : impl_(make_unique<Impl>())
    {
        // TODO check parameter but do not set amesos2 parameters
        set_parameters(params);
    }
    
    /**
     * Copy Constructor.
     * \param other an object Amesos2Solver
     */
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::Amesos2Solver(const Amesos2Solver &other)
    : impl_(make_unique<Impl>(*other.impl_)) {
        //FIXME
    }
    
    /**
     * Destructor.
     * \param na
     */
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::~Amesos2Solver() {}
    
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::Amesos2Solver() : impl_(make_unique<Impl>()) {}
    
    /**
     * update Method.
     * \param op a RCP pointer to Matrix
     * \param prec a RCP pointer to Matrix
     */
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op,
                                                         const std::shared_ptr<const Matrix> &prec)
    {
        PreconditionedSolver::update(op, prec);
        // set_problem(*op);
    }
    
    /**
     * update Method.
     * \param op a RCP pointer to Matrix
     * \param prec a RCP pointer to Matrix
     */
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op)
    {
        //preordering();
        // sym_factorization();
        // num_factorization();//TODO check output boolean       
        PreconditionedSolver::update(op);
        // set_problem(*op);
    }
    
    /**
     * apply Method - create the solver and and solve the system.
     * \param rhs, the right hand side vector
     * \param lhs, the left hand side vector
     * \return true
     */    
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::apply (const Vector &rhs, Vector &lhs) {
        using MatImplT = typename Matrix::Implementation::crs_mat_type;
        using VecImplT = typename Vector::Implementation::multi_vector_type;

        std::string solver_type = impl_->param_list_->sublist("UTOPIA", true).get("Solver Type", "LU");
        assert( Amesos2::query(solver_type) ); //check if the solver type specified in the xml is available
        
        impl_->solver_ = Amesos2::create<MatImplT, VecImplT>(
          solver_type.c_str(),
          raw_type(*this->get_operator()),
          raw_type(lhs),
          raw_type(rhs)
        );

        check_parameters();
        assert(!impl_->solver_.is_null());
        
        impl_->solver_->solve();
        return true;
    }
    
    
    /**
     * num_factorization Method - does the numeric factorization
     * \param na
     * \return bool
     */   
    // template <typename Matrix, typename Vector>
    // auto Amesos2Solver<Matrix, Vector, TRILINOS>::num_factorization() const 
    // { assert(!impl_->solver_.is_null());
    //     return impl_->solver_->numericFactorization(); }
    
    /**
     * sym_factorization Method - does the numeric factorization
     * \param na
     * \return bool
     */   
    // template <typename Matrix, typename Vector>
    // auto Amesos2Solver<Matrix, Vector, TRILINOS>::sym_factorization() const 
    // { assert(!impl_->solver_.is_null());
    //     return impl_->solver_->symbolicFactorization(); }
    
    /**
     * get_nnzLU Method.
     * \param na
     * \return int
     */   
    template <typename Matrix, typename Vector>
    int Amesos2Solver<Matrix, Vector, TRILINOS>::get_nnzLU() const { 
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().getNnzLU(); }
    
    /**
     * get_num_preorder Method - Returns the number of pre-orderings performed by the owning solver.
     * \param na
     * \return int
     */   
    template <typename Matrix, typename Vector>
    int     Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_preorder () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().getNumPreOrder(); }
    
    /**
     * get_num_sym_fact Method - Returns the number of symbolic factorizations performed by the owning solver.
     * \param na
     * \return int
     */   
    template <typename Matrix, typename Vector>
    int     Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_sym_fact () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().getNumSymbolicFact(); }
    
    /**
     * get_num_numeric_fact Method - Returns the number of numeric factorizations performed by the owning solver.
     * \param na
     * \return int
     */   
    template <typename Matrix, typename Vector>
    int     Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_numeric_fact () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().getNumNumericFact(); }
    
    /**
     * get_num_solve Method - Returns the number of solves performed by the owning solver.
     * \param na
     * \return int
     */   
    template <typename Matrix, typename Vector>
    int     Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_solve () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().getNumSolve(); } 
    /**
     * preordering_done Method - If true , then pre-ordering has been performed
     * \param na
     * \return int
     */   
    template <typename Matrix, typename Vector>
    bool    Amesos2Solver<Matrix, Vector, TRILINOS>::preordering_done () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().preOrderingDone(); }
    
    /**
     * sym_factorization_done Method - if true , then symbolic factorization has been performed.
     * \param na
     * \return int
     */   
    template <typename Matrix, typename Vector>
    bool    Amesos2Solver<Matrix, Vector, TRILINOS>::sym_factorization_done () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().symbolicFactorizationDone(); }
    
    /**
     * num_factorization_done Method -    If true , then numeric factorization has been performed.
     * \param na
     * \return int
     */   
    template <typename Matrix, typename Vector>
    bool    Amesos2Solver<Matrix, Vector, TRILINOS>::num_factorization_done () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().numericFactorizationDone(); }
    
    
    
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
    {
        set_preconditioner(); //(A) //FIXME
    }
    
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::set_preconditioner(const Matrix &precond)
    {//TODO useless for now
        bool direct_solver = impl_->param_list_->sublist("UTOPIA", true).template get<bool>("Direct Preconditioner", false);
        std::string dir_prec_type = impl_->param_list_->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");
        if ( direct_solver ) {
            #ifdef WITH_TRILINOS_IFPACK2
            impl_->ifpack_prec_ = Ifpack2::Factory::create<typename Impl::matrix_type>(dir_prec_type, precond.implementation().implementation_ptr());
            assert(!impl_->ifpack_prec_.is_null());
            impl_->ifpack_prec_->setParameters(impl_->param_list_->sublist(dir_prec_type, false));
            impl_->ifpack_prec_->initialize();
            impl_->ifpack_prec_->compute();
            //impl_->linear_problem->setLeftPrec(impl_->ifpack_prec_);
            #endif //WITH_TRILINOS_IFPACK2
        } else {
            #ifdef WITH_TRILINOS_MUELU
            // Multigrid Hierarchy
            impl_->muelu_prec_ = MueLu::CreateTpetraPreconditioner((
                                                                    Teuchos::RCP<typename Impl::OP>) precond.implementation().implementation_ptr(),
                                                                   impl_->param_list_->sublist("MueLu", false)
                                                                   );
            
            assert(!impl_->muelu_prec_.is_null());
            //impl_->linear_problem->setRightPrec(impl_->muelu_prec_);
            // #else
            assert(false);
            #endif //WITH_TRILINOS_MUELU
        }
    }
    
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::set_parameters(const Parameters params)
    {
        if(!params.param_file_name().empty()) {
            try {
                impl_->param_list_ = Teuchos::getParametersFromXmlFile(params.param_file_name());  //TODO this call should go in Param class for Trilinos together with the full param list
                auto& utopia_param_list = impl_->param_list_->sublist("UTOPIA", true);

                if( utopia_param_list.template get<bool>("Direct Solver", "true") )
                {
                 std::cout << "Using Direct Solvers " << std::endl;
                   if( Amesos2::query( utopia_param_list.get("Solver Type", "KLU2") ) ); //Amesos2::query returns true if the solver exists TODO print error
                   {
                    std::cout << "Solver Type: " << utopia_param_list.get("Solver Type", "KLU2")  << std::endl;
                   }
                }
                
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
    void Amesos2Solver<Matrix, Vector, TRILINOS>::check_parameters()
    {
           try {
Teuchos::RCP<Teuchos::ParameterList> tmp_param_list;
tmp_param_list.reset(new Teuchos::ParameterList(impl_->param_list_->sublist("Amesos2", true))); //TODO restrict to only Amesos2 input
impl_->solver_->setParameters(tmp_param_list);
            } catch(const std::exception &ex) {
                std::cerr << ex.what() << std::endl;
                assert(false);
                abort();
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
    
    // template <typename Matrix, typename Vector>
    // bool Amesos2Solver<Matrix, Vector, TRILINOS>::set_problem()
    // {
    //     impl_->linear_problem->setProblem();
    //     impl_->solver_ = impl_->amesos2_factory.create( impl_->param_list_->sublist("UTOPIA", true).get("Solver Type", "CG"), impl_->param_list_); //to change it to have the specialization
    //     impl_->solver_->setProblem(impl_->linear_problem);
    //     if (this->verbose()) { impl_->solver_->getCurrentParameters()->print(); }
    
    //     return true;
    // }
    
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::set_problem(Matrix &A)
    {
        // impl_->linear_problem->setProblem();
        set_preconditioner(); //TODO is it feasable??
        // impl_->solver_->setProblem(impl_->linear_problem);
        // if (this->verbose()) { impl_->solver_->getCurrentParameters()->print(); } //TODO print current parameters
        return true;
    }
    
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::set_preconditioner()//const std::shared_ptr<Preconditioner> &precond)
    {
        bool direct_solver = impl_->param_list_->sublist("UTOPIA", true).template get<bool>("Direct Preconditioner", false);
        std::string dir_prec_type = impl_->param_list_->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");
        
        if ( direct_solver )
        {
            // impl_->ifpack_prec_ = Ifpack2::Factory::create<Impl::matrix_type>(dir_prec_type, (Teuchos::RCP<Impl::OP>)precond->implementation().implementation_ptr()); //TODO
            
            #ifdef WITH_TRILINOS_IFPACK2
            assert(!impl_->ifpack_prec_.is_null());
            impl_->ifpack_prec_->setParameters(impl_->param_list_->sublist(dir_prec_type, false));
            impl_->ifpack_prec_->initialize();
            impl_->ifpack_prec_->compute();
            //            impl_->linear_problem->setLeftPrec(impl_->ifpack_prec_);
            #endif //WITH_TRILINOS_IFPACK2
            
        } else {
            #ifdef WITH_TRILINOS_MUELU
            // Multigrid Hierarchy
            //impl_->muelu_prec_ = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<Impl::OP>)precond->implementation().implementation_ptr(),    //TODO
            //                                                            impl_->param_list_->sublist("MueLu", false));
            assert(!impl_->muelu_prec_.is_null());
            //            impl_->linear_problem->setRightPrec(impl_->muelu_prec_);
            // #else
            assert(false);
            #endif //WITH_TRILINOS_MUELU
        }
    }
    
}  // namespace utopia

#endif //HAVE_AMESOS2_KOKKOS
#endif //UTOPIA_AMESOS2_IMPL_HPP

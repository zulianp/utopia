#ifndef UTOPIA_AMESOS2_IMPL_HPP
#define UTOPIA_AMESOS2_IMPL_HPP

#include "utopia_Amesos2_solver.hpp"

#include "utopia_make_unique.hpp"

#include <Amesos2_Factory.hpp>
#include <Amesos2_Solver.hpp>
#include "Amesos2.hpp"
#include "Amesos2_Meta.hpp"

//TODO remove from here
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>


// #ifdef HAVE_AMESOS2_TPETRA

//FIXME find right macros (these packages are optional in trilinos, they should be optional also in utopia)
// #define HAVE_AMESOS2_MUELU
// #define HAVE_AMESOS2_IFPACK2


// #ifdef HAVE_AMESOS2_MUELU
#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
// #else
// #warning "HAVE_AMESOS2_MUELU not defined"
// #endif //HAVE_AMESOS2_MUELU


// #ifdef HAVE_AMESOS2_IFPACK2
#include <Ifpack2_Factory.hpp>
// #else
// #warning "HAVE_AMESOS2_IFPACK2 not defined"
// #endif //HAVE_AMESOS2_IFPACK


namespace utopia {
/**
 * The Amesos2Solver class is the implementaion class
 *
 * \author Nur Aiman Fadel
 */
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

        typedef Tpetra::Vector<SC, LO, GO, NT> vec_type;
        typedef Tpetra::CrsMatrix<SC, LO, GO, NT> matrix_type;
        
        //typedef Amesos2::LinearProblem<SC, MV, OP> problem_type;
        typedef Amesos2::Solver<matrix_type, vec_type > solver_type;
// #ifdef HAVE_AMESOS2_IFPACK2
        typedef Ifpack2::Preconditioner<SC, LO, GO, NT> ifpack_prec_type;
// #endif //HAVE_AMESOS2_IFPACK
        
// #ifdef HAVE_AMESOS2_MUELU
        typedef MueLu::TpetraOperator<SC, LO, GO, NT> muelu_prec_type;
// #endif

// Members
        //Teuchos::RCP<problem_type> linear_problem;
        Teuchos::RCP<Teuchos::ParameterList> param_list_;
        //  auto& utopiaPL;// impl_->param_list_->sublist("UTOPIA", true);
        Teuchos::RCP<solver_type> solver_;
        
        //preconditioner
// #ifdef HAVE_AMESOS2_IFPACK2
        Teuchos::RCP<ifpack_prec_type> ifpack_prec_;
// #endif //HAVE_AMESOS2_IFPACK2
        
// #ifdef HAVE_AMESOS2_MUELU
        Teuchos::RCP<muelu_prec_type> muelu_prec_;
// #endif //HAVE_AMESOS2_MUELU
        
    };

  /**
   * Constructor that sets the solver's parameters.
   * \param params an object Parameters
   */
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::Amesos2Solver(Parameters params)
    : impl_(make_unique<Impl>())
    {
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

       std::string solver_type = impl_->param_list_->sublist("UTOPIA", true).get("Solver Type", "LU");
       assert( !Amesos2::query(solver_type) ); //check if the solver type specified in the xml is available

       impl_->solver_ = Amesos2::create<typename Impl::matrix_type, typename Impl::vec_type>
                        (solver_type, raw_type(*this->get_operator()), raw_type(lhs), raw_type(rhs));
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
// #ifdef HAVE_AMESOS2_IFPACK2
            impl_->ifpack_prec_ = Ifpack2::Factory::create<typename Impl::matrix_type>(dir_prec_type, precond.implementation().implementation_ptr());
            assert(!impl_->ifpack_prec_.is_null());
            impl_->ifpack_prec_->setParameters(impl_->param_list_->sublist(dir_prec_type, false));
            impl_->ifpack_prec_->initialize();
            impl_->ifpack_prec_->compute();
            //impl_->linear_problem->setLeftPrec(impl_->ifpack_prec_);
// #endif //HAVE_AMESOS2_IFPACK2
        } else {
// #ifdef HAVE_AMESOS2_MUELU
            // Multigrid Hierarchy
            impl_->muelu_prec_ = MueLu::CreateTpetraPreconditioner((
                                                                Teuchos::RCP<typename Impl::OP>) precond.implementation().implementation_ptr(),
                                                               impl_->param_list_->sublist("MueLu", false)
                                                               );
            
            assert(!impl_->muelu_prec_.is_null());
            //impl_->linear_problem->setRightPrec(impl_->muelu_prec_);
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
                impl_->param_list_ = Teuchos::getParametersFromXmlFile(params.param_file_name());
                impl_->solver_->setParameters(impl_->param_list_);
                assert( !Amesos2::query( impl_->param_list_->sublist("UTOPIA", true).get("Solver Type", "LU") ) );

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
            
// #ifdef HAVE_AMESOS2_IFPACK2
            assert(!impl_->ifpack_prec_.is_null());
            impl_->ifpack_prec_->setParameters(impl_->param_list_->sublist(dir_prec_type, false));
            impl_->ifpack_prec_->initialize();
            impl_->ifpack_prec_->compute();
//            impl_->linear_problem->setLeftPrec(impl_->ifpack_prec_);
// #endif //HAVE_AMESOS2_IFPACK2
            
        } else {
// #ifdef HAVE_AMESOS2_MUELU
            // Multigrid Hierarchy
            //impl_->muelu_prec_ = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<Impl::OP>)precond->implementation().implementation_ptr(),    //TODO
            //                                                            impl_->param_list_->sublist("MueLu", false));
            assert(!impl_->muelu_prec_.is_null());
//            impl_->linear_problem->setRightPrec(impl_->muelu_prec_);
// #else
            assert(false);
// #endif //HAVE_AMESOS2_MUELU
        }
    }
    
}  // namespace utopia

// #endif //HAVE_AMESOS2_TPETRA
#endif //UTOPIA_AMESOS2_IMPL_HPP

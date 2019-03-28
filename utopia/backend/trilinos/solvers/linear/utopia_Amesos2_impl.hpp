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
#include <Tpetra_Version.hpp>

//////



#ifdef HAVE_AMESOS2_KOKKOS

namespace utopia {
    /**
     * The Amesos2Solver class is the implementation class
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

#ifdef  KOKKOS_ENABLE_CUDA
        typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
        typedef cuda_node NT;
#elif defined  KOKKOS_ENABLE_ROCM //Kokkos::Compat::KokkosROCmWrapperNode doesn't exist
        typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::ROCm> rocm_node;
        typedef rocm_node NT;
#elif defined   KOKKOS_ENABLE_OPENMP
        typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
        typedef openmp_node NT;
#elif defined   KOKKOS_ENABLE_THREAD
        typedef Kokkos::Compat::KokkosThreadsWrapperNode thread_node;
#else
        typedef serial_node NT;
#endif

        typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
        typedef Tpetra::Operator<SC, LO, GO, NT> OP;

        typedef Tpetra::CrsMatrix<SC, LO, GO, NT> matrix_type;

        typedef Amesos2::Solver<matrix_type, multi_vector_type> solver_type;

        // Members
        Teuchos::RCP<Teuchos::ParameterList> amesos_list_;
        Teuchos::RCP<Teuchos::ParameterList> utopia_list_;

        Teuchos::RCP<solver_type> solver_;
        // Teuchos::RCP<matrix_type> matrix_;

        bool keep_symbolic_factorization;

        Impl()
        : keep_symbolic_factorization(false)
        {}
    };

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
     * update Method. - it replaces the matrix and does the preordering, symbolic and numericf actorization
     * \param op a RCP pointer to Matrix
     */
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op)
    {
        using MatImplT = typename Matrix::Implementation::crs_mat_type;
        using VecImplT = typename Vector::Implementation::multi_vector_type;

        DirectSolver<Matrix, Vector>::update(op);

        std::string solver_type = impl_->utopia_list_->get("Solver Type", "LU");
        assert( Amesos2::query(solver_type) ); //check if the solver type specified in the xml is available

        //CHECK IF SOLVER HAS CHANGED

        bool first = false;
        if(impl_->solver_.is_null()) {
            // impl_->matrix_ = raw_type(*op);
            impl_->solver_ = Amesos2::create<MatImplT, VecImplT>(solver_type,
                                                                 raw_type(*op)
                                                                 );
            first = true;
        }

        // possible options are CLEAN, PREORDERING, SYMBFACT, NUMFACT, SOLVE
        impl_->solver_->setA(
            raw_type(*op),
            (
              impl_->keep_symbolic_factorization ?
              Amesos2::SYMBFACT :
              Amesos2::CLEAN
              )
            );

        // with SYMBFACT you keep the symbolic factorization

        if(!impl_->keep_symbolic_factorization || first) {
            preordering();
            sym_factorization();
        }


        num_factorization();
    }

    /**
     * apply Method - create the solver and and solve the system.
     * \param rhs, the right hand side vector
     * \param lhs, the left hand side vector
     * \return true
     */
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::apply(const Vector &rhs, Vector &lhs) {
        assert(!impl_->solver_.is_null());

        if(impl_->solver_.is_null()) {
            utopia_error("solver not initialized. Call update first");
            return false;
        }

        impl_->solver_->setX(raw_type(lhs));
        impl_->solver_->setB(raw_type(rhs));

        check_parameters();
        // impl_->matrix_.reset();
        impl_->solver_->solve();//TODO does preordering, sym_factorization, num_factorization ?
        return true;
    }

    /**
     * get_preordering_done Method - If true , then pre-ordering has been performed
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool    Amesos2Solver<Matrix, Vector, TRILINOS>::get_preordering_done () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().preOrderingDone(); }

    /**
     * get_sym_factorization_done Method - if true , then symbolic factorization has been performed.
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool    Amesos2Solver<Matrix, Vector, TRILINOS>::get_sym_factorization_done () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().symbolicFactorizationDone(); }

    /**
     * get_num_factorization_done Method -    If true , then numeric factorization has been performed.
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool    Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_factorization_done () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().numericFactorizationDone(); }

    /**
     * preordering Method - does the preordering of the matrix entries
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::preordering()
    { assert(!impl_->solver_.is_null());
        impl_->solver_->preOrdering();
        return get_preordering_done ();
    }

    /**
     * num_factorization Method - does the numeric factorization
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::num_factorization()
    { assert(!impl_->solver_.is_null());
        impl_->solver_->numericFactorization();
        return get_num_factorization_done ();
    }

    /**
     * sym_factorization Method - does the numeric factorization
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::sym_factorization()
    { assert(!impl_->solver_.is_null());
        impl_->solver_->symbolicFactorization();
        return get_sym_factorization_done ();
    }

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
    int Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_solve () const {
        assert(!impl_->solver_.is_null());
        return impl_->solver_->getStatus().getNumSolve();
    }

    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::read(Input &in)
    {
        DirectSolver<Matrix, Vector>::read(in);

        in.get("keep-symbolic-factorization", impl_->keep_symbolic_factorization);
    }

    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::print_usage(std::ostream &os) const
    {

    }

    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::read_xml(const std::string &path)
    {
        if(!path.empty()) {
            try {
                Teuchos::RCP<Teuchos::ParameterList> tmp_param_list;
                tmp_param_list = Teuchos::getParametersFromXmlFile(path);  //TODO this call should go in Param class for Trilinos together with the full param list

                impl_->amesos_list_.reset(new Teuchos::ParameterList(tmp_param_list->sublist("Amesos2", true)));
                impl_->utopia_list_.reset(new Teuchos::ParameterList(tmp_param_list->sublist("UTOPIA", true)));
                //impl_->utopia_list_->template get<bool>("Direct Solver", "true")
                // Amesos2::query( impl_->utopia_list_->get("Solver Type", "KLU2") ) // Amesos2::query returns true if the solver exists TODO print error

                //TODO: move validation to param class

                /*        Teuchos::RCP<Teuchos::ParameterList> correct_list;

                 Teuchos::setStringToIntegralParameter<Amesos2::EPhase>("Direct Sol. Phase to use", "CLEAN",
                 "Update options for Matrix A",
                 Teuchos::tuple<std::string>("CLEAN","PREORDERING","SYMBFACT","NUMFACT","SOLVE"),
                 Teuchos::tuple<std::string>("Start from scratch",
                 "Keep Preordering",
                 "Keep Symbolic Factorization",
                 "Keep Numeric Factorization",
                 "Keep Solution"),
                 Teuchos::tuple<Amesos2::EPhase>(Amesos2::CLEAN,
                 Amesos2::PREORDERING,
                 Amesos2::SYMBFACT,
                 Amesos2::NUMFACT,
                 Amesos2::SOLVE),
                 correct_list.getRawPtr());


                 Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = correct_list->getEntry("Direct Sol. Phase to use").validator();
                 impl_->utopia_list_->getEntry("Direct Sol. Phase to use").setValidator(validator);*/

            } catch(const std::exception &ex) {
                std::cerr << ex.what() << std::endl;
                assert(false);
                Utopia::Abort();
            }
        } else {
            //TODO use default paramlist
        }
    }

    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::check_parameters(){
        try {
            impl_->solver_->setParameters(impl_->amesos_list_);
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

}  // namespace utopia

#endif //HAVE_AMESOS2_KOKKOS
#endif //UTOPIA_AMESOS2_IMPL_HPP

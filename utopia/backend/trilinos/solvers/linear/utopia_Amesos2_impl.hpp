#ifndef UTOPIA_AMESOS2_IMPL_HPP
#define UTOPIA_AMESOS2_IMPL_HPP

#include "utopia_Amesos2_solver.hpp"
#include "utopia_Base.hpp"
#include "utopia_make_unique.hpp"

#ifdef UTOPIA_WITH_TRILINOS_AMESOS2
#include <Amesos2.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

namespace utopia {
    /**
     * The Amesos2Solver class is the implementation class
     *
     * \author Nur Aiman Fadel
     */
    template <typename Matrix, typename Vector>
    class Amesos2Solver<Matrix, Vector, TRILINOS>::Impl {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using LocalSizeType = typename Traits<Vector>::LocalSizeType;
        using Node = typename Traits<Vector>::Node;

        using MultiVectorType = typename Vector::MultiVectorType;
        using CrsMatrixType = typename Matrix::CrsMatrixType;
        using SolverType = Amesos2::Solver<CrsMatrixType, MultiVectorType>;

        // Members
        Teuchos::RCP<Teuchos::ParameterList> param_list;
        Teuchos::RCP<SolverType> amesos_solver;

        Impl() : param_list(Teuchos::parameterList()) {}
    };

    /**
     * Copy Constructor.
     * \param other an object Amesos2Solver
     */
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::Amesos2Solver(const Amesos2Solver &other)
        : solver_type_(other.solver_type_), impl_(make_unique<Impl>(*other.impl_)) {
        // FIXME
    }

    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::Amesos2Solver(const std::string &type)
        : solver_type_(type), impl_(make_unique<Impl>()) {
        // check if the solver type specified as constructor parameter is available
        assert(Amesos2::query(solver_type_));

        // retrieve solver config from xml file
        const auto tmp_param_list =
            Teuchos::getParametersFromXmlFile(Utopia::instance().get("data_path") + "/xml/UTOPIA_amesos.xml");
        impl_->param_list.reset(new Teuchos::ParameterList(tmp_param_list->sublist(solver_type_, false)));
    }

    /**
     * Destructor.
     * \param na
     */
    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS>::~Amesos2Solver() {}

    /**
     * update Method. - it replaces the matrix and does the preordering, symbolic and numericf actorization
     * \param op a RCP pointer to Matrix
     */
    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op) {
        using CrsMatrixType = typename Matrix::CrsMatrixType;
        using MultiVectorType = typename Vector::MultiVectorType;

        DirectSolver<Matrix, Vector>::update(op);

        // CHECK IF SOLVER HAS CHANGED

        bool first = false;
        if (impl_->amesos_solver.is_null()) {
            // impl_->matrix_ = raw_type(*op);
            impl_->amesos_solver = Amesos2::create<CrsMatrixType, MultiVectorType>(solver_type_, raw_type(*op));
            first = true;
        }

        // possible options are CLEAN, PREORDERING, SYMBFACT, NUMFACT, SOLVE
        impl_->amesos_solver->setA(raw_type(*op), (keep_symbolic_factorization ? Amesos2::SYMBFACT : Amesos2::CLEAN));

        // with SYMBFACT you keep the symbolic factorization

        if (!keep_symbolic_factorization || first) {
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
        assert(!impl_->amesos_solver.is_null());

        if (impl_->amesos_solver.is_null()) {
            utopia_error("solver not initialized. Call update first");
            return false;
        }

        impl_->amesos_solver->setX(raw_type(lhs));
        impl_->amesos_solver->setB(raw_type(rhs));

        try {
            impl_->amesos_solver->setParameters(impl_->param_list);
        } catch (const std::exception &ex) {
            std::cerr << ex.what() << std::endl;
            assert(0);
        }
        impl_->amesos_solver->solve();  // TODO does preordering, sym_factorization, num_factorization ?
        return true;
    }

    /**
     * get_preordering_done Method - If true , then pre-ordering has been performed
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::get_preordering_done() const {
        assert(!impl_->amesos_solver.is_null());
        return impl_->amesos_solver->getStatus().preOrderingDone();
    }

    /**
     * get_sym_factorization_done Method - if true , then symbolic factorization has been performed.
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::get_sym_factorization_done() const {
        assert(!impl_->amesos_solver.is_null());
        return impl_->amesos_solver->getStatus().symbolicFactorizationDone();
    }

    /**
     * get_num_factorization_done Method -    If true , then numeric factorization has been performed.
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_factorization_done() const {
        assert(!impl_->amesos_solver.is_null());
        return impl_->amesos_solver->getStatus().numericFactorizationDone();
    }

    /**
     * preordering Method - does the preordering of the matrix entries
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::preordering() {
        assert(!impl_->amesos_solver.is_null());
        impl_->amesos_solver->preOrdering();
        return get_preordering_done();
    }

    /**
     * num_factorization Method - does the numeric factorization
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::num_factorization() {
        assert(!impl_->amesos_solver.is_null());
        impl_->amesos_solver->numericFactorization();
        return get_num_factorization_done();
    }

    /**
     * sym_factorization Method - does the numeric factorization
     * \param na
     * \return bool
     */
    template <typename Matrix, typename Vector>
    bool Amesos2Solver<Matrix, Vector, TRILINOS>::sym_factorization() {
        assert(!impl_->amesos_solver.is_null());
        impl_->amesos_solver->symbolicFactorization();
        return get_sym_factorization_done();
    }

    /**
     * get_nnzLU Method.
     * \param na
     * \return int
     */
    template <typename Matrix, typename Vector>
    int Amesos2Solver<Matrix, Vector, TRILINOS>::get_nnzLU() const {
        assert(!impl_->amesos_solver.is_null());
        return impl_->amesos_solver->getStatus().getNnzLU();
    }

    /**
     * get_num_preorder Method - Returns the number of pre-orderings performed by the owning solver.
     * \param na
     * \return int
     */
    template <typename Matrix, typename Vector>
    int Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_preorder() const {
        assert(!impl_->amesos_solver.is_null());
        return impl_->amesos_solver->getStatus().getNumPreOrder();
    }

    /**
     * get_num_sym_fact Method - Returns the number of symbolic factorizations performed by the owning solver.
     * \param na
     * \return int
     */
    template <typename Matrix, typename Vector>
    int Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_sym_fact() const {
        assert(!impl_->amesos_solver.is_null());
        return impl_->amesos_solver->getStatus().getNumSymbolicFact();
    }

    /**
     * get_num_numeric_fact Method - Returns the number of numeric factorizations performed by the owning solver.
     * \param na
     * \return int
     */
    template <typename Matrix, typename Vector>
    int Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_numeric_fact() const {
        assert(!impl_->amesos_solver.is_null());
        return impl_->amesos_solver->getStatus().getNumNumericFact();
    }

    /**
     * get_num_solve Method - Returns the number of solves performed by the owning solver.
     * \param na
     * \return int
     */
    template <typename Matrix, typename Vector>
    int Amesos2Solver<Matrix, Vector, TRILINOS>::get_num_solve() const {
        assert(!impl_->amesos_solver.is_null());
        return impl_->amesos_solver->getStatus().getNumSolve();
    }

    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::read(Input &in) {
        DirectSolver<Matrix, Vector>::read(in);
        in.get("keep_symbolic_factorization", keep_symbolic_factorization);
    }

    template <typename Matrix, typename Vector>
    void Amesos2Solver<Matrix, Vector, TRILINOS>::print_usage(std::ostream &os) const {
        DirectSolver<Matrix, Vector>::print_usage(os);
        this->print_param_usage(
            os, "keep_symbolic_factorization", "bool", "flag to keep the symbolic factorization", "false");
    }

    template <typename Matrix, typename Vector>
    Amesos2Solver<Matrix, Vector, TRILINOS> *Amesos2Solver<Matrix, Vector, TRILINOS>::clone() const {
        return new Amesos2Solver(*this);
    }

}  // namespace utopia

#endif  // UTOPIA_AMESOS2_IMPL_HPP
#endif  // UTOPIA_WITH_TRILINOS_AMESOS2
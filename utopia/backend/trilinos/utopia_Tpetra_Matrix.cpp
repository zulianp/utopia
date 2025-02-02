#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_Matrix_impl.hpp"

#include "utopia_Allocations.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_Instance.hpp"
#include "utopia_Logger.hpp"
// #include "utopia_kokkos_ParallelEach.hpp"
//  #include "utopia_trilinos_Each_impl.hpp"

#include "utopia_trilinos_Utils.hpp"

#include <MatrixMarket_Tpetra.hpp>
#include <TpetraExt_MatrixMatrix_def.hpp>
#include <Tpetra_RowMatrixTransposer_decl.hpp>

#include <array>
#include <iterator>

// FIXME
// - crs matrix has problematic behaviour when adding to off-procs entries once assembled in finite element assembly
// routines

namespace utopia {

    void TpetraMatrix::set(const SizeType &row, const SizeType &col, const Scalar &value) {
        m_utopia_status_once(
            "> TpetraMatrix::set does what is supposed to do with respect to the edsl. "
            "However it might be slow because of the double trilinos call and branching.");

        if (implementation().replaceGlobalValues(row, 1, &value, &col) != 1) {
            implementation().insertGlobalValues(row, 1, &value, &col);
        }
    }

    void TpetraMatrix::add(const SizeType &row, const SizeType &col, const Scalar &value) {
        m_utopia_status_once(
            "> TpetraMatrix::add does what is supposed to do with respect to the edsl. "
            "However it might be slow because of the double trilinos call and branching.");

        if (implementation().sumIntoGlobalValues(row, 1, &value, &col, false) != 1) {
            implementation().insertGlobalValues(row, 1, &value, &col);
        }
    }

    void TpetraMatrix::c_set(const SizeType &row, const SizeType &col, const Scalar &value) {
        m_utopia_status_once(
            "> TpetraMatrix::set does what is supposed to do with respect to the edsl. "
            "However it might be slow because of the double trilinos call and branching.");

        if (implementation().replaceGlobalValues(row, 1, &value, &col) != 1) {
            implementation().insertGlobalValues(row, 1, &value, &col);
        }
    }

    void TpetraMatrix::c_add(const SizeType &row, const SizeType &col, const Scalar &value) {
        m_utopia_status_once(
            "> TpetraMatrix::add does what is supposed to do with respect to the edsl. "
            "However it might be slow because of the double trilinos call and branching.");

        if (implementation().sumIntoGlobalValues(row, 1, &value, &col, false) != 1) {
            implementation().insertGlobalValues(row, 1, &value, &col);
        }
    }

    // FIXME make faster version by storing view?
    TpetraMatrix::Scalar TpetraMatrix::get(const SizeType &row, const SizeType &col) const {
#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
        CrsMatrixType::local_inds_host_view_type cols;
        CrsMatrixType::values_host_view_type values;
#else
        Teuchos::ArrayView<const LocalSizeType> cols;
        Teuchos::ArrayView<const Scalar> values;
#endif

        assert(implementation().isLocallyIndexed());

        auto local_col = col - implementation().getDomainMap()->getMinGlobalIndex();

        auto rr = row_range();
        implementation().getLocalRowView(row - rr.begin(), cols, values);

        auto it = std::lower_bound(cols.data(), cols.data() + cols.size(), local_col);

        if (it == cols.data() + cols.size()) {
            return 0.;
        }

        // assert(it != std::end(cols));

        std::size_t index = std::distance(cols.data(), it);

        assert(cols[index] == local_col);

        return values[index];
    }

    void TpetraMatrix::multiply(const TpetraVector &vec, TpetraVector &result) const {
        assert(vec.size() == cols());

        if (result.is_null()) {
            result.init(mat_->getRowMap());
            // result.owner_ = true;
        } else if (!result.implementation().getMap()->isSameAs(*mat_->getRowMap())) {
            result.init(mat_->getRowMap());
            // result.owner_ = true;
        }
        try {
            mat_->apply(vec.implementation(), result.implementation());
        } catch (const std::exception &ex) {
            utopia::out() << ex.what() << std::endl;
            assert(false);
        }

        assert(result.size() == rows());
    }

    void TpetraMatrix::transpose_multiply(const TpetraVector &vec, TpetraVector &result) const {
        assert(mat_->hasTransposeApply());

        if (result.is_null()) {
            result.init(mat_->getDomainMap());
            // result.owner_ = true;
        } else if (!result.implementation().getMap()->isSameAs(*mat_->getDomainMap())) {
            result.init(mat_->getDomainMap());
            // result.owner_ = true;
        }
        try {
            mat_->apply(vec.implementation(), result.implementation(), Teuchos::TRANS);
        } catch (const std::exception &ex) {
            utopia::out() << ex.what() << std::endl;
            assert(false);
        }
    }

    void TpetraMatrix::multiply(const TpetraMatrix &right, TpetraMatrix &result) const {
        multiply(false, false, right, result);
    }

    void TpetraMatrix::multiply(const Scalar &alpha, const TpetraMatrix &B, TpetraMatrix &C) const {
        multiply(B, C);
        C.scale(alpha);
    }

    void TpetraMatrix::transpose_multiply(const TpetraMatrix &right, TpetraMatrix &result) const {
        multiply(true, false, right, result);
    }

    // result op(*this) * op
    void TpetraMatrix::multiply(const bool transpose_this,
                                const bool transpose_right,
                                const TpetraMatrix &right,
                                TpetraMatrix &result) const {
        m_utopia_status_once("TpetraMatrix::multiply Proper thing to do would be to check if the maps are compatible");
        // IMPROVEME
        result.mat_.reset();

        assert(!transpose_right);

        assert(transpose_this ||
               (this->local_size().get(1) == right.local_size().get(0) && this->size().get(1) == right.size().get(0)));
        assert(!transpose_this ||
               (this->local_size().get(0) == right.local_size().get(0) && this->size().get(0) == right.size().get(0)));

        if (transpose_this) {
            assert(!right.implementation().getDomainMap().is_null());

            UTOPIA_REPORT_ALLOC("TpetraMatrix::multiply");
            result.mat_ = Teuchos::rcp(new CrsMatrixType(implementation().getDomainMap(),
                                                         // col_map,
                                                         0
#if (TRILINOS_MAJOR_MINOR_VERSION < 130100 || !UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
                                                         ,
                                                         Tpetra::StaticProfile
#endif
                                                         ));

        } else {
            UTOPIA_REPORT_ALLOC("TpetraMatrix::multiply");
            result.mat_ = Teuchos::rcp(new CrsMatrixType(implementation().getRowMap(),
                                                         // col_map,
                                                         0
#if (TRILINOS_MAJOR_MINOR_VERSION < 130100 || !UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
                                                         ,
                                                         Tpetra::StaticProfile
#endif
                                                         ));
        }

        result.owner_ = true;

        try {
            Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
            params->set("Optimize Storage", true);

            Tpetra::MatrixMatrix::Multiply(this->implementation(),
                                           transpose_this,
                                           right.implementation(),
                                           transpose_right,
                                           result.implementation(),
                                           true,
                                           "TpetraMatrix::multiply",
                                           params);

        } catch (const std::exception &ex) {
            utopia::out() << ex.what() << std::endl;
            assert(false);
        }
    }

    bool TpetraMatrix::is_assembled() const { return raw_type()->isFillComplete(); }

    void TpetraMatrix::transpose(TpetraMatrix &mat) const {
        try {
            Tpetra::RowMatrixTransposer<Scalar, LocalSizeType, SizeType, Node> transposer(mat_);

            auto temp = transposer.createTranspose();

            mat.mat_ = temp;
            mat.owner_ = true;
            assert(is_valid(true));
        } catch (const std::exception &ex) {
            utopia::out() << ex.what() << std::endl;
            assert(false);
        }
    }

    void TpetraMatrix::axpy(const Scalar &alpha, const TpetraMatrix &x) {
        try {
            mat_ = Tpetra::MatrixMatrix::add(alpha, false, x.implementation(), 1., false, implementation());

            owner_ = true;
        } catch (const std::exception &ex) {
            utopia::out() << ex.what() << std::endl;
            assert(false);
        }
    }

    void TpetraMatrix::finalize() {
        try {
            if (init_) {
                assert(!init_->domain_map.is_null());
                assert(!init_->range_map.is_null());

                implementation().fillComplete(init_->domain_map, init_->range_map);
            } else {
                implementation().fillComplete(implementation().getDomainMap(), implementation().getRangeMap());
            }

        } catch (const std::exception &ex) {
            utopia::out() << ex.what() << std::endl;
            assert(false);
        }
    }

    void TpetraMatrix::crs_init(const RCPCommType &comm,
                                std::size_t rows_local,
                                std::size_t cols_local,
                                Tpetra::global_size_t rows_global,
                                Tpetra::global_size_t cols_global,
                                std::size_t nnz_x_row) {
        // Trilinos has more distribution options than petsc and it does not require to have
        // a column operator structure as petsc

        RCPMapType row_map;
        const int index_base = 0;

        if (rows_local == static_cast<std::size_t>(INVALID_INDEX)) {
            // NEW SIZE
            const SizeType rows_local = utopia::decompose(comm_, rows_global);
            row_map.reset(new MapType(rows_global, rows_local, index_base, comm));
        } else {
            row_map.reset(new MapType(rows_global, rows_local, index_base, comm));
        }

        if (cols_global == Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()) {
            Tpetra::global_size_t send_buff = cols_local;
            cols_global = 0;
            Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &send_buff, &cols_global);
        }

        // auto col_map = Teuchos::rcp(new MapType(cols_global, index_base, comm, Tpetra::LocallyReplicated));
        // mat_.reset(new CrsMatrixType(row_map, col_map, nnz_x_row, Tpetra::StaticProfile));

        UTOPIA_REPORT_ALLOC("TpetraMatrix::crs_init");
        mat_.reset(new CrsMatrixType(row_map,
                                     nnz_x_row
#if (TRILINOS_MAJOR_MINOR_VERSION < 130100 || !UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
                                     ,
                                     Tpetra::StaticProfile
#endif

                                     ));
        owner_ = true;

        init_ = std::make_shared<InitStructs>();
        if (cols_local == static_cast<std::size_t>(INVALID_INDEX)) {
            // NEW SIZE
            const SizeType cols_local = utopia::decompose(comm_, cols_global);
            init_->domain_map.reset(new MapType(cols_global, cols_local, index_base, comm));
        } else {
            init_->domain_map.reset(new MapType(cols_global, cols_local, index_base, comm));
        }

        init_->range_map = row_map;
    }

    void TpetraMatrix::crs_init(const RCPCommType &comm,
                                std::size_t rows_local,
                                std::size_t cols_local,
                                Tpetra::global_size_t rows_global,
                                Tpetra::global_size_t cols_global,
                                const Teuchos::ArrayRCP<size_t> &rowPtr,
                                const Teuchos::ArrayRCP<LocalSizeType> &cols,
                                const Teuchos::ArrayRCP<Scalar> &values) {
        RCPMapType row_map;
        RCPMapType col_map;
        const int index_base = 0;
        if (rows_local == static_cast<std::size_t>(INVALID_INDEX)) {
            assert(cols_local == static_cast<std::size_t>(INVALID_INDEX));

            const SizeType rows_local_auto = utopia::decompose(comm_, rows_global);

            row_map.reset(new MapType(rows_global, rows_local_auto, index_base, comm));
            col_map.reset(new MapType(cols_global, index_base, comm));
            //~ Kokkos::View<LO*> colInds("Column Map", cols_global);
            //~ Kokkos::parallel_for(cols_global, KOKKOS_LAMBDA(size_t i) { colInds(i) = i; });
            //~ col_map.reset(new MapType(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
            // Kokkos::Compat::getConstArrayView(colInds), index_base, comm));
        } else {
            // see for a distributed example
            // https://github.com/trilinos/Trilinos/blob/master/packages/tpetra/core/example/Lesson07-Kokkos-Fill/04_tpetra.cpp
            // assert(false && "Sparse distributed matrix assembly with CRS structures is not implemented yet.");

            row_map.reset(new MapType(rows_global, rows_local, index_base, comm));
            col_map.reset(new MapType(cols_global, index_base, comm));
        }

        try {
            UTOPIA_REPORT_ALLOC("TpetraMatrix::crs_init");
            mat_.reset(new CrsMatrixType(row_map, col_map, rowPtr, cols, values));
            owner_ = true;

            init_ = std::make_shared<InitStructs>();
            if (cols_local == static_cast<std::size_t>(INVALID_INDEX)) {
                const SizeType cols_local_auto = utopia::decompose(comm_, cols_global);
                init_->domain_map.reset(new MapType(cols_global, cols_local_auto, index_base, comm));
            } else {
                init_->domain_map.reset(new MapType(cols_global, cols_local, index_base, comm));
            }
            init_->range_map = row_map;

            finalize();
        } catch (const std::exception &ex) {
            utopia::out() << ex.what() << std::endl;
            assert(false);
            throw ex;
        }
    }

    void TpetraMatrix::crs_identity(const RCPCommType &comm,
                                    std::size_t rows_local,
                                    std::size_t cols_local,
                                    Tpetra::global_size_t rows_global,
                                    Tpetra::global_size_t cols_global,
                                    const Scalar factor) {
        crs_init(comm, rows_local, cols_local, rows_global, cols_global, 1.);

        write_lock();

        Range r = row_range();
        const std::size_t r_begin = r.begin();
        const std::size_t r_end = r.end();

        auto cols = init_->domain_map->getGlobalNumElements();

        for (auto i = r_begin; i < r_end; ++i) {
            if (i < cols) {
                set(i, i, factor);
            } else {
                break;
            }
        }

        write_unlock();
    }

    void TpetraMatrix::build_diag(TpetraVector &d) const {
        const bool is_row_min = this->size().get(0) <= this->size().get(1);
        SizeType n = (is_row_min) ? this->size().get(0) : this->size().get(1);

        if (d.is_null() || d.size() != n) {
            m_utopia_warning_once("TpetraMatrix::get_diag Assuming row <= col");

            if (is_row_min) {
                d.init(implementation().getRowMap());
            } else {
                d.init(implementation().getDomainMap());
            }
        }

        implementation().getLocalDiagCopy(d.implementation());
    }

    void TpetraMatrix::diag(const TpetraMatrix &mat) {
        TpetraVector d(mat.comm());
        mat.build_diag(d);
        diag(d);
    }

    void TpetraMatrix::diag(const TpetraVector &d) {
        auto ls = d.local_size();
        auto gs = d.size();

        crs_init(d.communicator(), ls, ls, gs, gs, 1);

        auto r = d.range();
        auto data = d.implementation().getData();

        assert(!data.is_null());

        write_lock();

        LocalSizeType index = 0;

        for (auto i = r.begin(); i < r.end(); ++i) {
            set(i, i, data[index++]);
        }

        write_unlock();
    }

    bool TpetraMatrix::read(const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const std::string &path) {
        std::ifstream is;
        is.open(path.c_str());

        if (!is.good()) {
            return false;
        }

        try {
            // https://people.sc.fsu.edu/~jburkardt/data/mm/mm.html
            mat_ = Tpetra::MatrixMarket::Reader<CrsMatrixType>::readSparse(is, comm);
        } catch (std::exception &ex) {
            is.close();
            utopia::out() << ex.what() << std::endl;
            return false;
        }

        is.close();
        return !mat_.is_null();
    }

    bool TpetraMatrix::write(const std::string &path) const {
        if (mat_.is_null()) {
            return false;
        }

        try {
            Tpetra::MatrixMarket::Writer<CrsMatrixType>::writeSparseFile(path, mat_, "mat", "", false);
        } catch (const std::exception &ex) {
            utopia::out() << ex.what() << std::endl;
            return false;
        }

        return true;
    }

    bool TpetraMatrix::is_valid(const bool verbose) const {
        if (mat_.is_null()) {
            if (verbose) {
                std::cerr << "is_null" << std::endl;
            }
            return false;
        }

        auto comm = communicator();

        if (comm->getSize() == 1) {
            if (local_size() != size()) {
                if (verbose) {
                    std::cerr << "local_size() != size()" << std::endl;
                    std::cerr << local_size() << " != " << size() << std::endl;
                    std::cerr << "this indicates inconsistent domain_map with respect to the col_map" << std::endl;
                }

                return false;
            }
        }

        return true;
    }

    TpetraMatrix::Scalar TpetraMatrix::norm2() const { return implementation().getFrobeniusNorm(); }

    TpetraMatrix::Scalar TpetraMatrix::sum() const {
        TpetraVector vec(this->comm()), row_sum(this->comm());
        vec.values(layout(this->communicator(), this->local_size().get(1), this->size().get(1)), 1.);
        this->multiply(vec, row_sum);
        return row_sum.sum();
    }

    TpetraMatrix::Scalar TpetraMatrix::norm_infty() const {
        Scalar ret = 0;
        return parallel_reduce_values(AbsMax(), Max(), ret);
    }

    TpetraMatrix::Scalar TpetraMatrix::norm1() const { return parallel_reduce_values(AbsPlus(), Plus(), 0); }

    TpetraMatrix::SizeType TpetraMatrix::rows() const {
        if (is_null()) {
            return 0;
        }

        if (implementation().isFillComplete()) {
            return implementation().getGlobalNumRows();
        }
        assert(!implementation().getRowMap().is_null());
        return implementation().getRowMap()->getGlobalNumElements();
    }

    TpetraMatrix::SizeType TpetraMatrix::cols() const {
        if (is_null()) {
            return 0;
        }

        if (implementation().isFillComplete()) {
            return implementation().getGlobalNumCols();
        }
        if (implementation().getDomainMap().is_null()) {
            assert(!init_->domain_map.is_null());
            return init_->domain_map->getGlobalNumElements();
        }
        return implementation().getDomainMap()->getGlobalNumElements();
    }

    TpetraMatrix::SizeType TpetraMatrix::local_rows() const {
        if (is_null()) {
            return 0;
        }

        assert(!implementation().getRowMap().is_null());

#if UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE == 1
        return implementation().getRowMap()->getLocalNumElements();
#else
        return implementation().getRowMap()->getNodeNumElements();
#endif
    }

    TpetraMatrix::SizeType TpetraMatrix::local_cols() const {
        if (is_null()) {
            return 0;
        }

        if (implementation().getDomainMap().is_null()) {
            assert(!init_->domain_map.is_null());
#if UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE == 1
            return init_->domain_map->getLocalNumElements();
#else
            return init_->domain_map->getNodeNumElements();
#endif
        }

#if UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE == 1
        return implementation().getDomainMap()->getLocalNumElements();
#else
        return implementation().getDomainMap()->getNodeNumElements();
#endif
    }

    void TpetraMatrix::clear() {
        mat_.reset();
        owner_ = true;
        init_.reset();
    }

    void TpetraMatrix::set_zero_rows(const IndexSet &index, const Scalar &diag) {
        auto &impl = implementation();
        auto rr = row_range();

        auto col_map = impl.getColMap()->getLocalMap();
        auto row_map = impl.getRowMap()->getLocalMap();

#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
        auto local_mat = raw_type()->getLocalMatrixHost();
#else
        auto local_mat = raw_type()->getLocalMatrix();
#endif

        for (auto i_global : index) {
            if (!rr.inside(i_global)) {
                std::cerr << "[Error] index out of range " << i_global << " not in " << rr << std::endl;
                assert(rr.inside(i_global));
                continue;
            }

            auto i = i_global - rr.begin();
            auto row = local_mat.row(i);
            auto n_values = row.length;

            for (decltype(n_values) k = 0; k < n_values; ++k) {
                auto &val = row.value(k);
                const auto col = row.colidx(k);

                if (row_map.getGlobalElement(i) == col_map.getGlobalElement(col)) {
                    val = diag;
                } else {
                    val = 0.;
                }
            }
        }
    }

    void TpetraMatrix::select(const IndexSet & /*row_index*/,
                              const IndexSet & /*col_index*/,
                              TpetraMatrix & /*result*/) const {
        assert(false && "IMPLEMENT ME");
    }

    void TpetraMatrix::transform(const Sqrt &op) { aux_transform(op); }

    void TpetraMatrix::transform(const Pow2 &op) { aux_transform(op); }

    void TpetraMatrix::transform(const Log &op) { aux_transform(op); }

    void TpetraMatrix::transform(const Exp &op) { aux_transform(op); }

    void TpetraMatrix::transform(const Cos &op) { aux_transform(op); }

    void TpetraMatrix::transform(const Sin &op) { aux_transform(op); }

    void TpetraMatrix::transform(const Abs &op) { aux_transform(op); }

    void TpetraMatrix::transform(const Minus &op) { aux_transform(op); }

    void TpetraMatrix::transform(const Pow &op) { aux_transform(op); }

    void TpetraMatrix::transform(const Reciprocal<Scalar> &op) { aux_transform(op); }

    void TpetraMatrix::swap(TpetraMatrix &x) {
        using std::swap;
        swap(comm_, x.comm_);
        swap(mat_, x.mat_);
        swap(owner_, x.owner_);
    }

    void TpetraMatrix::scale(const Scalar &alpha) {
        // Why???
        write_lock();
        implementation().scale(alpha);
        write_unlock();
    }

    TpetraMatrix::Scalar TpetraMatrix::dot(const TpetraMatrix &other) const {
        if (raw_type() == other.raw_type()) {
            const Scalar ret = norm2();
            return ret * ret;
        }

        assert(false && "IMPLEMENT ME");
        return -1.0;
    }

    //(y := alpha * A * x + beta * y)
    void TpetraMatrix::gemv(const bool transpose,
                            const Scalar &alpha,
                            const TpetraVector &x,
                            const Scalar &beta,
                            TpetraVector &y) const {
        if (alpha == 0.0) {
            y.scale(beta);
            return;
        }

        if (beta == 0.0) {
            if (transpose) {
                transpose_multiply(x, y);
            } else {
                multiply(x, y);
            }

            y.scale(alpha);
            return;
        }

        // BAD we have a temporary here
        TpetraVector temp(this->comm());

        if (transpose) {
            transpose_multiply(x, temp);
        } else {
            multiply(x, temp);
        }

        temp.scale(alpha);
        temp.axpy(beta, y);
        y.copy(temp);
    }

    bool TpetraMatrix::equals(const TpetraMatrix &other, const Scalar &tol) const {
        if (raw_type() == other.raw_type()) {
            return true;
        }

        if (rows() != other.rows() || cols() != other.cols()) {
            return false;
        }

        TpetraMatrix diff = *this;
        diff.axpy(-1.0, other);
        return diff.norm_infty() <= tol;
    }

    void TpetraMatrix::build_from_structure(const TpetraMatrix &rhs) {
        UTOPIA_REPORT_ALLOC("TpetraMatrix::build_from_structure");
        const auto &rhs_ptr = rhs.raw_type();
        owner_ = true;
        mat_.reset(new CrsMatrixType(rhs_ptr->getCrsGraph()));
    }

    void TpetraMatrix::shift_diag(const TpetraVector &d) {
        if (empty()) {
            diag(d);
            return;
        }

        auto d_view = const_view_device(d);

        // FIXME?
        this->transform_ijv(UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &v)->Scalar {
            return (i == j) ? (v + d_view.get(i)) : v;
        });
    }

    void TpetraMatrix::set_diag(const TpetraVector &d) {
        auto d_view = const_view_device(d);

        // FIXME?
        this->transform_ijv(UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &v)->Scalar {
            return (i == j) ? (d_view.get(i)) : v;
        });
    }

    void TpetraMatrix::diag_scale_left(const TpetraVector & /*d*/) { assert(false && "IMPLEMENT ME"); }

}  // namespace utopia

#include "utopia_petsc_CrsView.hpp"

#include "utopia_Reporter.hpp"

#include "utopia_petsc_Matrix.hpp"

namespace utopia {

    class PetscCrsView::Impl {
    public:
        Impl() = default;
        Impl(Mat raw_mat) { init(raw_mat); }
        ~Impl() { destroy(); }

        void set(Mat raw_mat) {
            destroy();
            init(raw_mat);
        }

        Mat raw_mat{nullptr};
        PetscInt rows{0};
        PetscInt cols{0};
        const PetscInt *ia{nullptr};
        const PetscInt *ja{nullptr};
        PetscScalar *array{nullptr};
        PetscBool done{PETSC_FALSE};
        PetscErrorCode err{0};

        void destroy() {
            if (raw_mat) {
                if (PetscMatrix::is_block(raw_mat)) {
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 13, 4)
                    MatSeqBAIJRestoreArray(raw_mat, &array);
#else
                    utopia::err() << "[Error] MATBAIJ not supported\n";
#endif
                } else {
                    MatSeqAIJRestoreArray(raw_mat, &array);
                }

                err = MatRestoreRowIJ(raw_mat, 0, PETSC_FALSE, PETSC_FALSE, &rows, &ia, &ja, &done);
                assert(err == 0);
                UTOPIA_UNUSED(err);
                raw_mat = nullptr;
                array = nullptr;
            }
        }

        void init(Mat raw_mat) {
            this->raw_mat = raw_mat;

            err = MatGetRowIJ(raw_mat, 0, PETSC_FALSE, PETSC_FALSE, &rows, &ia, &ja, &done);

            assert(err == 0);
            assert(done == PETSC_TRUE);

            err = MatGetSize(raw_mat, nullptr, &cols);
            assert(err == 0);

            if (!done) {
                utopia::err()
                    << "PetscMatrix::read_petsc_seqaij_impl(const Op &op): MatGetRowIJ failed to provide what "
                       "was asked.\n";
                abort();
            }

            if (PetscMatrix::is_block(raw_mat)) {
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 13, 4)
                MatSeqBAIJGetArray(raw_mat, &array);
#else
                utopia::err() << "[Error] MATBAIJ not supported\n";
#endif
            } else {
                MatSeqAIJGetArray(raw_mat, &array);
            }
        }
    };

    void PetscCrsView::clear() { impl_->destroy(); }

    PetscCrsView::PetscCrsView(Mat raw_mat) : impl_(std::make_shared<Impl>(raw_mat)) {}
    PetscCrsView::PetscCrsView() : impl_(std::make_shared<Impl>()) {}
    PetscCrsView::~PetscCrsView() = default;

    void PetscCrsView::set(Mat raw_mat) { impl_->set(raw_mat); }

    PetscInt PetscCrsView::nnz() const { return impl_->ia[impl_->rows]; }

    PetscInt PetscCrsView::cols() const { return impl_->cols; }

    PetscInt PetscCrsView::rows() const { return impl_->rows; }

    PetscCrsView::RowView PetscCrsView::row(const PetscInt row) const {
        if (empty()) {
            return RowView();
        }

        PetscInt row_begin = impl_->ia[row];
        PetscInt n = impl_->ia[row + 1] - row_begin;
        return RowView(n, &impl_->ja[row_begin], &impl_->array[row_begin]);
    }

    ArrayView<const PetscInt> PetscCrsView::row_ptr() const { return ArrayView<const PetscInt>(impl_->ia, rows() + 1); }

    ArrayView<const PetscInt> PetscCrsView::colidx() const { return ArrayView<const PetscInt>(impl_->ja, nnz()); }

    ArrayView<const PetscScalar> PetscCrsView::values() const {
        return ArrayView<const PetscScalar>(impl_->array, nnz());
    }

    ArrayView<PetscScalar> PetscCrsView::values() { return ArrayView<PetscScalar>(impl_->array, nnz()); }

    void PetscCrsView::describe(std::ostream &os) const {
        os << rows() << " " << cols() << '\n';
        os << "nnz: " << nnz() << '\n';
    }

    bool PetscCrsView::empty() const { return !impl_ || !impl_->raw_mat; }

    PetscCrsView crs_view(PetscMatrix &mat) {
        assert(!mat.is_mpi());
        return mat.raw_type();
    }

    PetscCrsView crs_view(const PetscMatrix &mat) {
        assert(!mat.is_mpi());
        return mat.raw_type();
    }

    void views_host(PetscMatrix &mat, PetscCrsView &d, PetscCrsView &o) {
        if (mat.is_mpi()) {
            Mat d_, o_;
            MatMPIAIJGetSeqAIJ(mat.raw_type(), &d_, &o_, nullptr);
            d.set(d_);
            o.set(o_);
        } else {
            d.set(mat.raw_type());
            o.clear();
        }
    }

}  // namespace utopia

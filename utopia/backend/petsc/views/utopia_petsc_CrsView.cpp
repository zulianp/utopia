#include "utopia_petsc_CrsView.hpp"

#include "utopia_Reporter.hpp"

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

    private:
        void destroy() {
            if (raw_mat) {
                MatSeqAIJRestoreArray(raw_mat, &array);
                err = MatRestoreRowIJ(raw_mat, 0, PETSC_FALSE, PETSC_FALSE, &rows, &ia, &ja, &done);
                assert(err == 0);
                UTOPIA_UNUSED(err);
                raw_mat = nullptr;
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

            MatSeqAIJGetArray(raw_mat, &array);
        }
    };

    PetscCrsView::PetscCrsView(Mat raw_mat) : impl_(std::make_shared<Impl>(raw_mat)) {}
    PetscCrsView::PetscCrsView() : impl_(std::make_shared<Impl>()) {}
    PetscCrsView::~PetscCrsView() = default;

    void PetscCrsView::set(Mat raw_mat) { impl_->set(raw_mat); }

    PetscInt PetscCrsView::nnz() const { return impl_->ia[impl_->rows]; }

    PetscInt PetscCrsView::cols() const { return impl_->rows; }

    PetscInt PetscCrsView::rows() const { return impl_->rows; }

    ArrayView<const PetscInt> PetscCrsView::row_ptr() const { return ArrayView<const PetscInt>(impl_->ia, rows() + 1); }

    ArrayView<const PetscInt> PetscCrsView::colidx() const { return ArrayView<const PetscInt>(impl_->ja, nnz()); }

    ArrayView<const PetscScalar> PetscCrsView::values() const {
        return ArrayView<const PetscScalar>(impl_->array, nnz());
    }

    ArrayView<PetscScalar> PetscCrsView::values() { return ArrayView<PetscScalar>(impl_->array, nnz()); }

}  // namespace utopia

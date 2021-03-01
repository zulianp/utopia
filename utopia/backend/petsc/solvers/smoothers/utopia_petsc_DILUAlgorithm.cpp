#include "utopia_petsc_DILUAlgorithm.hpp"

#include "utopia_make_unique.hpp"

#include "utopia_petsc.hpp"

#include "utopia_petsc_CrsView.hpp"
#include "utopia_petsc_ILUDecompose.hpp"

#include "utopia_DILUDecompose_impl.hpp"

namespace utopia {

    class DILUAlgorithm<PetscMatrix, PetscVector>::Impl {
    public:
        using ScalarView = utopia::ArrayView<PetscScalar>;
        using IndexView = utopia::ArrayView<const PetscInt>;

        bool update(const PetscMatrix &mat) {
            if (mat.comm().size() == 1) {
                local_mat.copy(mat);
            } else {
                local_block_view(mat, local_mat);
            }

            crs_view.set(local_mat.raw_type());
            return dilu.update(crs_view);
        }

        PetscMatrix local_mat;
        PetscCrsView crs_view;
        DILUAlgorithm<PetscCrsView, ScalarView> dilu;
    };

    void DILUAlgorithm<PetscMatrix, PetscVector>::init() {
        impl_ = utopia::make_unique<Impl>();
        if (!impl_) {
            utopia::err() << "[Error] DILUAlgorithm<PetscMatrix, PetscVector>::construct: impl_ == nullptr\n";
            Utopia::Abort();
        } else {
            // utopia::out() << "[Status] init called\n";
        }
    }

    DILUAlgorithm<PetscMatrix, PetscVector>::DILUAlgorithm() : Super() { init(); }

    DILUAlgorithm<PetscMatrix, PetscVector>::~DILUAlgorithm() {}

    bool DILUAlgorithm<PetscMatrix, PetscVector>::update(const PetscMatrix &mat) {
        assert(impl_);
        if (impl_) {
            return impl_->update(mat);
        } else {
            utopia::err() << "[Error] DILUAlgorithm<PetscMatrix, PetscVector>::update(...): impl_ == nullptr\n";
            Utopia::Abort();
            return false;
        }
    }

    void DILUAlgorithm<PetscMatrix, PetscVector>::apply(const PetscVector &b, PetscVector &x) {
        auto x_view = local_view_device(x);
        auto b_view = local_view_device(const_cast<PetscVector &>(b));

        auto x_array = x_view.array();

        assert(impl_);

        if (impl_) {
            impl_->dilu.apply(b_view.array(), x_array);
        } else {
            utopia::err() << "[Error] DILUAlgorithm<PetscMatrix, PetscVector>::apply(...): impl_ == nullptr\n";
            Utopia::Abort();
        }
    }

    void DILUAlgorithm<PetscMatrix, PetscVector>::read(Input &) {}

}  // namespace utopia

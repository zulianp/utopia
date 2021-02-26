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
            local_block_view(mat, local_mat);
            crs_view.set(local_mat.raw_type());
            return dilu.update(crs_view);
        }

        PetscMatrix local_mat;
        PetscCrsView crs_view;
        DILUAlgorithm<PetscCrsView, ScalarView> dilu;
    };

    DILUAlgorithm<PetscMatrix, PetscVector>::DILUAlgorithm() : impl_(utopia::make_unique<Impl>()) {}

    DILUAlgorithm<PetscMatrix, PetscVector>::~DILUAlgorithm() {}

    bool DILUAlgorithm<PetscMatrix, PetscVector>::update(const PetscMatrix &mat) { return impl_->update(mat); }

    void DILUAlgorithm<PetscMatrix, PetscVector>::apply(const PetscVector &b, PetscVector &x) {
        auto x_view = local_view_device(x);
        auto b_view = local_view_device(const_cast<PetscVector &>(b));

        auto x_array = x_view.array();

        impl_->dilu.apply(b_view.array(), x_array);
    }

    void DILUAlgorithm<PetscMatrix, PetscVector>::read(Input &) {}

}  // namespace utopia

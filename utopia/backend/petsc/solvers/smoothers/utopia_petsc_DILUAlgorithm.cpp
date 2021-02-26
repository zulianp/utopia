#include "utopia_petsc_DILUAlgorithm.hpp"

#include "utopia_make_unique.hpp"

#include "utopia_petsc.hpp"

namespace utopia {

    class DILUAlgorithm<PetscMatrix, PetscVector>::Impl {
    public:
        using ScalarView = utopia::ArrayView<PetscScalar>;
        using IndexView = utopia::ArrayView<const PetscInt>;

        PetscMatrix local_mat;
        CRSMatrix<ScalarView, IndexView> crs_mat;
    };

    DILUAlgorithm<PetscMatrix, PetscVector>::DILUAlgorithm() : impl_(utopia::make_unique<Impl>()) {}

    DILUAlgorithm<PetscMatrix, PetscVector>::~DILUAlgorithm() {}

    bool DILUAlgorithm<PetscMatrix, PetscVector>::update(const PetscMatrix &mat) {
        local_block_view(mat, impl_->local_mat);

        return false;
    }

    void DILUAlgorithm<PetscMatrix, PetscVector>::apply(const PetscVector &b, PetscVector &x) {}

    void DILUAlgorithm<PetscMatrix, PetscVector>::read(Input &) {}

}  // namespace utopia
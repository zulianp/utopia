
#include "utopia_petsc_Eval_Chop.hpp"
#include "utopia_petsc_Matrix_impl.hpp"

namespace utopia {

    void Chop<PetscMatrix, PETSC>::apply(PetscMatrix &mat, const Scalar &eps) {
        // each_transform(mat, [eps](const SizeType &, const SizeType &, const Scalar &v) -> Scalar {
        // 	return std::abs(v) < eps ? 0.0 : v;
        // });

        mat.transform_values([eps](const Scalar &v) -> Scalar { return std::abs(v) < eps ? 0.0 : v; });
    }

    void ChopSmallerThan<PetscMatrix, PETSC>::apply(PetscMatrix &mat, const Scalar &eps) {
        // each_transform(mat, [eps](const SizeType &, const SizeType &, const Scalar &v) -> Scalar {
        // 	return v < eps ? 0.0 : v;
        // });

        mat.transform_values([eps](const Scalar &v) -> Scalar { return v < eps ? 0.0 : v; });
    }

    void ChopGreaterThan<PetscMatrix, PETSC>::apply(PetscMatrix &mat, const Scalar &eps) {
        // each_transform(mat, [eps](const SizeType &, const SizeType &, const Scalar &v) -> Scalar {
        // 	return v > eps ? 0.0 : v;
        // });

        mat.transform_values([eps](const Scalar &v) -> Scalar { return v > eps ? 0.0 : v; });
    }

}  // namespace utopia

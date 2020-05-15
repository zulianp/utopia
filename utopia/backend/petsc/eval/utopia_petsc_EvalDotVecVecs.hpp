#ifndef UTOPIA_EVAL_PETSC_DOT_VEC_VECS_HPP
#define UTOPIA_EVAL_PETSC_DOT_VEC_VECS_HPP

#include <vector>
#include "utopia_DotVecVecs.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
    template <class Vector>
    class EvalDots<Vector, PETSC> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        static void apply(const Vector &v1,
                          const std::vector<std::shared_ptr<Vector> > &vectors,
                          std::vector<Scalar> &results);

        static void apply(const Vector &v1, const std::vector<Vector> &vectors, std::vector<Scalar> &results);

        static void apply(const Vector &v11,
                          const Vector &v12,
                          Scalar &result1,
                          const Vector &v21,
                          const Vector &v22,
                          Scalar &result2);

        static void apply(const Vector &v11,
                          const Vector &v12,
                          Scalar &result1,
                          const Vector &v21,
                          const Vector &v22,
                          Scalar &result2,
                          const Vector &v31,
                          const Vector &v32,
                          Scalar &result3);
    };

    template <class Vector>
    class EvalNorm2s<Vector, PETSC> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        static void apply(const Vector &v1, const Vector &v2, Scalar &result1, Scalar &result2);

        static void apply(const Vector &v1,
                          const Vector &v2,
                          const Vector &v3,
                          Scalar &result1,
                          Scalar &result2,
                          Scalar &result3);
    };

}  // namespace utopia

#endif  // UTOPIA_EVAL_PETSC_DOT_VEC_VECS_HPP

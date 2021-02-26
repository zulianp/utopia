#ifndef UTOPIA_DILU_DECOMPOSE_HPP
#define UTOPIA_DILU_DECOMPOSE_HPP

#include "utopia_CRSMatrix.hpp"
#include "utopia_ILUDecompose.hpp"
#include "utopia_Traits.hpp"

#include "utopia_CrsMatrixIndexer.hpp"

#include <vector>

namespace utopia {

    template <class Matrix, class Vector>
    class DILUAlgorithm final : public ILUAlgorithm<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;

        bool decompose(const Matrix &in, std::vector<Scalar> &d);

        bool update(const Matrix &mat) override;
        void apply(const Vector &b, Vector &x) override;
        void read(Input &) override;

    private:
        std::shared_ptr<const Matrix> mat_;
        std::vector<Scalar> d_, L_inv_b_;

        CrsDiagIndexer<SizeType> diag_idx_;
        CrsTransposeIndexer<SizeType> transpose_idx_;
    };

}  // namespace utopia

#endif  // UTOPIA_DILU_DECOMPOSE_HPP

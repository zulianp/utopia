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

        bool update(const Matrix &mat) override;
        void apply(const Vector &b, Vector &x) override;
        void read(Input &) override;
        inline DILUAlgorithm *clone() const override { return new DILUAlgorithm(); }

    private:
        std::shared_ptr<const Matrix> mat_;
        std::vector<Scalar> d_, L_inv_b_;

        CrsDiagIndexer<SizeType> diag_idx_;
        CrsTransposeIndexer<SizeType> transpose_idx_;
    };

    template <class Matrix, class Vector, int BlockSize>
    class BlockDILUAlgorithm final : public ILUAlgorithm<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;
        using ArrayViewT = utopia::ArrayView<const Scalar, DYNAMIC_SIZE, DYNAMIC_SIZE>;
        using BlockView = utopia::TensorView<ArrayViewT, 2>;
        using VectorView = utopia::TensorView<ArrayView<Scalar, DYNAMIC_SIZE>, 1>;
        using ConstVectorView = utopia::TensorView<ArrayView<const Scalar, DYNAMIC_SIZE>, 1>;

        static const int BlockSize_2 = BlockSize * BlockSize;
        using Block = utopia::StaticMatrix<Scalar, BlockSize, BlockSize>;

        bool update(const Matrix &mat) override;
        void apply(const Vector &b, Vector &x) override;
        void read(Input &) override;

        inline BlockDILUAlgorithm *clone() const override { return new BlockDILUAlgorithm(); }

    private:
        std::shared_ptr<const Matrix> mat_;
        std::vector<Block> d_;
        std::vector<Scalar> L_inv_b_;

        CrsDiagIndexer<SizeType> diag_idx_;
        CrsTransposeIndexer<SizeType> transpose_idx_;
    };

}  // namespace utopia

#endif  // UTOPIA_DILU_DECOMPOSE_HPP

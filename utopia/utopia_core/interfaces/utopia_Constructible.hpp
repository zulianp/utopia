#ifndef UTOPIA_CONSTRUCTIBLE_HPP
#define UTOPIA_CONSTRUCTIBLE_HPP

#include "utopia_Size.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    template <class Tensor, int Order_ = Traits<Tensor>::Order>
    class Constructible {};

    template <class Tensor>
    class SparseConstructible {
    public:
        using Traits_ = utopia::Traits<Tensor>;
        using Scalar = typename Traits_::Scalar;
        using SizeType = typename Traits_::SizeType;
        using LocalSizeType = typename Traits_::LocalSizeType;
        using MatrixLayout = typename Traits_::MatrixLayout;

        virtual ~SparseConstructible() = default;

        virtual void sparse(const MatrixLayout &layout, const SizeType &nnz_d_block, const SizeType &nnz_o_block) = 0;
        virtual void identity(const MatrixLayout &layout, const Scalar &diag = 1.0) = 0;
    };

    template <class Tensor>
    class DenseConstructible {
    public:
        using Traits_ = utopia::Traits<Tensor>;
        using Scalar = typename Traits_::Scalar;
        using SizeType = typename Traits_::SizeType;
        using LocalSizeType = typename Traits_::LocalSizeType;
        using MatrixLayout = typename Traits_::MatrixLayout;

        virtual ~DenseConstructible() = default;

        virtual void dense(const MatrixLayout &layout, const Scalar &val = 0.0) = 0;
        virtual void dense_identity(const MatrixLayout &layout, const Scalar &diag = 1.0) = 0;
    };

    template <class Tensor>
    class Constructible<Tensor, 2> : public SparseConstructible<Tensor>, public DenseConstructible<Tensor> {
    public:
        using Traits_ = utopia::Traits<Tensor>;
        using Scalar = typename Traits_::Scalar;
        using SizeType = typename Traits_::SizeType;
        using LocalSizeType = typename Traits_::LocalSizeType;
        using MatrixLayout = typename Traits_::MatrixLayout;

        ~Constructible() override = default;
    };

    template <class Tensor>
    class Constructible<Tensor, 1> {
    public:
        using Traits_ = utopia::Traits<Tensor>;
        using Scalar = typename Traits_::Scalar;
        using SizeType = typename Traits_::SizeType;
        using LocalSizeType = typename Traits_::LocalSizeType;
        using Layout = typename Traits_::Layout;

        virtual ~Constructible() = default;
        virtual void values(const Layout &l, const Scalar &value) = 0;
        virtual void zeros(const Layout &l) = 0;
    };
}  // namespace utopia

#endif  // UTOPIA_CONSTRUCTIBLE_HPP

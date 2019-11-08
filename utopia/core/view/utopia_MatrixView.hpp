#ifndef UTOPIA_KOKKOS_MATRIX_VIEW_HPP
#define UTOPIA_KOKKOS_MATRIX_VIEW_HPP

#include "utopia_Tensor.hpp"
#include "utopia_Algorithms.hpp"
#include <utility>

namespace utopia {

    template<class ArrayView2D_>
    class MatrixView final : public Tensor<MatrixView<ArrayView2D_>, 2> {
    public:
        using ArrayView2D = ArrayView2D_;
        using Scalar   = typename Traits<ArrayView2D>::Scalar;
        using SizeType = typename Traits<ArrayView2D>::SizeType;
        
        using Super = utopia::Tensor<MatrixView, 2>;
        using Super::Super;

        inline std::string get_class() const override
        {
            return "MatrixView";
        }

        template<class... Args>
        UTOPIA_FUNCTION MatrixView(Args && ...args)
        : view_(std::forward<Args>(args)...)
        {}

        template<class Expr>
        UTOPIA_FUNCTION MatrixView(const Expression<Expr> &expr)
        {
            //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::construct_eval(expr.derived());
        }

        template<class Expr>
        UTOPIA_INLINE_FUNCTION MatrixView &operator=(const Expression<Expr> &expr)
        {
            //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::assign_eval(expr.derived());
            return *this;
        }

        template<class OtherArrayView2D>
        UTOPIA_FUNCTION void copy(const MatrixView<OtherArrayView2D> &other)
        {
            UTOPIA_DEVICE_ASSERT(rows() == other.rows());
            UTOPIA_DEVICE_ASSERT(cols() == other.cols());

            device::copy(other.view_, view_);
        }

        UTOPIA_FUNCTION void assign(const MatrixView &other) override
        {
            copy(other);
        }

        UTOPIA_FUNCTION void assign(MatrixView &&other) override
        {
            view_ = std::move(other.view_);
        }

        UTOPIA_INLINE_FUNCTION ArrayView2D &raw_type() { return view_; }
        UTOPIA_INLINE_FUNCTION const ArrayView2D &raw_type() const { return view_; }

        UTOPIA_INLINE_FUNCTION SizeType rows() const
        {
            return device::extent(view_, 0);
        }

        UTOPIA_INLINE_FUNCTION SizeType cols() const
        {
            return device::extent(view_, 1);
        }

        UTOPIA_INLINE_FUNCTION const Scalar &get(const SizeType &i, const SizeType &j) const
        {
            return view_(i, j);
        }

        UTOPIA_INLINE_FUNCTION void set(const SizeType &i, const SizeType &j, const Scalar &value)
        {
            view_(i, j) = value;
        }

        UTOPIA_INLINE_FUNCTION void add(const SizeType &i, const SizeType &j, const Scalar &value)
        {
            view_(i, j) += value;
        }

        UTOPIA_INLINE_FUNCTION void scale(const Scalar &alpha)
        {
            device::scale(alpha, view_);
        }

        UTOPIA_INLINE_FUNCTION void axpy(const Scalar &alpha, const MatrixView &x)
        {
            device::axpy(alpha, x.view_, view_);
        }

        UTOPIA_INLINE_FUNCTION Scalar dot(const MatrixView &other) const
        {
            return device::dot(view_, other.view_);
        }

        template<class RightView, class ResultView>
        UTOPIA_INLINE_FUNCTION void multiply(
            const VectorView<RightView> &right,
            VectorView<ResultView> &result) const
        {
            UTOPIA_DEVICE_ASSERT(!result.is_alias(right));
            UTOPIA_DEVICE_ASSERT(result.size() == rows());
            UTOPIA_DEVICE_ASSERT(right.size() == cols());

            device::mv(view_, right.raw_type(), result.raw_type());
        }

        template<class RightView, class ResultView>
        UTOPIA_INLINE_FUNCTION void multiply(const MatrixView<RightView> &right, MatrixView<ResultView> &result) const
        {
            UTOPIA_DEVICE_ASSERT(!result.is_alias(right));
            UTOPIA_DEVICE_ASSERT(!is_alias(result));
            UTOPIA_DEVICE_ASSERT(rows() == result.rows());
            UTOPIA_DEVICE_ASSERT(right.cols() == result.cols());
            UTOPIA_DEVICE_ASSERT(cols() == right.rows());

            device::mm(view_, right.view_, result.view_);
        }

        UTOPIA_INLINE_FUNCTION void set(const Scalar &alpha)
        {
            device::fill(alpha, view_);
        }

        UTOPIA_FUNCTION MatrixView(const ArrayView2D &view)
        : view_(view) {}

        inline void describe() const
        {
            const SizeType r = rows();
            const SizeType c = cols();

            for(SizeType i = 0; i < r; ++i) {
                for(SizeType j = 0; j < c; ++j) {
                    std::cout << get(i, j) << "\t";
                }

                std::cout << "\n";
            }
        }

        UTOPIA_INLINE_FUNCTION void wrap(const ArrayView2D &view)
        {
            view_ = view;
        }

        UTOPIA_INLINE_FUNCTION bool is_alias(const MatrixView &other) const
        {
            return &(view_(0, 0)) == &(other.view_(0, 0));
        }

        
        template<class OtherView>
        UTOPIA_INLINE_FUNCTION constexpr static bool is_alias(const MatrixView<OtherView> &)
        {
            return false;
        }

    private:
        ArrayView2D view_;

        UTOPIA_FUNCTION MatrixView(const MatrixView &) {
            UTOPIA_DEVICE_ASSERT(false);
        }
    };

}

#endif //UTOPIA_KOKKOS_MATRIX_VIEW_HPP

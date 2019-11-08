#ifndef UTOPIA_KOKKOS_VECTOR_VIEW_HPP
#define UTOPIA_KOKKOS_VECTOR_VIEW_HPP

#include "utopia_Tensor.hpp"
#include "utopia_Algorithms.hpp"

namespace utopia {

    template<class ArrayView_>
    class VectorView final : public Tensor<VectorView<ArrayView_>, 1> {
    public:
        using ArrayView = ArrayView_;
        using Scalar   = typename Traits<ArrayView>::Scalar;
        using SizeType = typename Traits<ArrayView>::SizeType;

        using Super = utopia::Tensor<VectorView, 1>;
        using Super::Super;

        inline std::string get_class() const override
        {
            return "VectorView";
        }

        
        template<class... Args>
        UTOPIA_FUNCTION VectorView(Args && ...args)
        : view_(std::forward<Args>(args)...)
        {}


        template<class Expr>
        UTOPIA_FUNCTION VectorView(const Expression<Expr> &expr)
        {
            //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::construct_eval(expr.derived());
        }

        template<class Expr>
        UTOPIA_INLINE_FUNCTION VectorView &operator=(const Expression<Expr> &expr)
        {
            //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::assign_eval(expr.derived());
            return *this;
        }

        template<class OtherArrayView>
        UTOPIA_FUNCTION void copy(const VectorView<OtherArrayView> &other)
        {
            UTOPIA_DEVICE_ASSERT(size() == other.size());

            // const SizeType n = size();

            // for(SizeType i = 0; i < n; ++i) {
            //     view_[i] = other.get(i);
            // }

            device::copy(other.view_, view_);
        }

        template<class OtherArrayView>
        UTOPIA_FUNCTION void assign(const VectorView<OtherArrayView> &other)
        {
            copy(other);
        }

        UTOPIA_FUNCTION void assign(const VectorView &other) override
        {
            copy(other);
        }

        UTOPIA_FUNCTION void assign(VectorView &&other) override
        {
            view_ = std::move(other.view_);
        }

        UTOPIA_INLINE_FUNCTION ArrayView &raw_type() { return view_; }
        UTOPIA_INLINE_FUNCTION const ArrayView &raw_type() const { return view_; }

        UTOPIA_INLINE_FUNCTION SizeType size() const
        {
            return view_.size();
        }

        UTOPIA_INLINE_FUNCTION const Scalar &get(const SizeType &i) const
        {
            return view_[i];
        }

        UTOPIA_INLINE_FUNCTION void set(const SizeType &i, const Scalar &value)
        {
            view_[i] = value;
        }

        UTOPIA_INLINE_FUNCTION void add(const SizeType &i, const Scalar &value)
        {
            view_[i] += value;
        }

        UTOPIA_INLINE_FUNCTION void scale(const Scalar &alpha)
        {
            device::scale(alpha, view_);
        }

        template<class OtherArrayView>
        UTOPIA_INLINE_FUNCTION void axpy(const Scalar &alpha, const VectorView<OtherArrayView> &x)
        {
            return device::axpy(alpha, x.view_, view_);
        }

        template<class OtherArrayView>
        UTOPIA_INLINE_FUNCTION Scalar dot(const VectorView<OtherArrayView> &other) const
        {
            return device::dot(view_, other.view_);
        }

        UTOPIA_INLINE_FUNCTION void set(const Scalar &alpha)
        {
            device::fill(alpha, view_);
        }

        UTOPIA_FUNCTION VectorView(const ArrayView &view)
        : view_(view) {}

        inline void describe() const
        {
            const SizeType n = size();
            for(SizeType i = 0; i < n; ++i) {
                std::cout << get(i) << std::endl;
            }
        }

        UTOPIA_INLINE_FUNCTION void wrap(const ArrayView &view)
        {
            view_ = view;
        }

        UTOPIA_INLINE_FUNCTION bool is_alias(const VectorView &other) const
        {
            return &(view_[0]) == &(other.view_[0]);
        }

        template<class OtherView>
        UTOPIA_INLINE_FUNCTION constexpr static bool is_alias(const VectorView<OtherView> &)
        {
            return false;
        }

        template<class OtherView>
        inline bool equals(const VectorView<OtherView> &other, const Scalar &tol) const
        {
            if(size() != other.size()) return false;
            return device::approxeq(view_, other.view_, tol);
        }

    private:
        ArrayView view_;

        UTOPIA_FUNCTION VectorView(const VectorView &other) : view_(other.view_) {
            UTOPIA_DEVICE_ASSERT(false);
        }
    };

}

#endif //UTOPIA_KOKKOS_VECTOR_VIEW_HPP

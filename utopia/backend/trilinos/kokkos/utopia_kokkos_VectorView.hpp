#ifndef UTOPIA_KOKKOS_VECTOR_VIEW_HPP
#define UTOPIA_KOKKOS_VECTOR_VIEW_HPP

#include "utopia_Tensor.hpp"
#include "utopia_kokkos_Base.hpp"

#include <Kokkos_View.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_scal.hpp>

namespace utopia {

    template<class KokkosView_>
    class VectorView final : public Tensor<VectorView<KokkosView_>, 1> {
    public:
        using KokkosView = KokkosView_;
        using Scalar = typename KokkosView::value_type;
        using SizeType = std::size_t;

        using Super = utopia::Tensor<VectorView, 1>;
        using Super::Super;

        inline std::string get_class() const override
        {
            return "VectorView";
        }

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

        UTOPIA_FUNCTION void assign(const VectorView &other) override
        {
            UTOPIA_DEVICE_ASSERT(size() == other.size());

            Kokkos::parallel_for("VectorView::assign", view_.size(), KOKKOS_LAMBDA (const int& i) {
                view_(i) = other.view_(i);
            });
        }

        UTOPIA_FUNCTION void assign(VectorView &&other) override
        {
            view_ = std::move(other.view_);
        }

        UTOPIA_INLINE_FUNCTION KokkosView &raw_type() { return view_; }
        UTOPIA_INLINE_FUNCTION const KokkosView &raw_type() const { return view_; }

        UTOPIA_INLINE_FUNCTION void resize(const SizeType &n)
        {
            Kokkos::resize(view_, n);
        }

        UTOPIA_INLINE_FUNCTION SizeType size() const
        {
            return view_.size();
        }

        UTOPIA_INLINE_FUNCTION const Scalar &get(const SizeType &i) const
        {
            return view_(i);
        }

        UTOPIA_INLINE_FUNCTION void set(const SizeType &i, const Scalar &value)
        {
            view_(i) = value;
        }

        UTOPIA_INLINE_FUNCTION void add(const SizeType &i, const Scalar &value)
        {
            view_(i) += value;
        }

        UTOPIA_INLINE_FUNCTION void scale(const Scalar &alpha)
        {
            KokkosBlas::scal(view_,alpha,view_);
        }

        UTOPIA_INLINE_FUNCTION void axpy(const Scalar &alpha, const VectorView &x)
        {
            KokkosBlas::axpy(alpha,x.view_,view_);
        }

        UTOPIA_INLINE_FUNCTION Scalar dot(const VectorView &other) const
        {
            return KokkosBlas::dot(view_, other.view_);
        }

        UTOPIA_INLINE_FUNCTION void set(const Scalar &alpha)
        {
            KokkosBlas::fill(view_, alpha);
        }

        UTOPIA_FUNCTION VectorView(const KokkosView &view)
        : view_(view) {}

        inline void describe() const
        {
            const SizeType n = size();
            for(SizeType i = 0; i < n; ++i) {
                std::cout << get(i) << std::endl;
            }
        }

        UTOPIA_INLINE_FUNCTION void wrap(const KokkosView &view)
        {
            view_ = view;
        }

        UTOPIA_INLINE_FUNCTION bool is_alias(const VectorView &other) const
        {
            return &(view_(0)) == &(other.view_(0));
        }

    private:
        KokkosView view_;

        UTOPIA_FUNCTION VectorView(const VectorView &) {
            UTOPIA_DEVICE_ASSERT(false);
        }
    };

}

#endif //UTOPIA_KOKKOS_VECTOR_VIEW_HPP

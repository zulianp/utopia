#ifndef UTOPIA_KOKKOS_MATRIX_VIEW_HPP
#define UTOPIA_KOKKOS_MATRIX_VIEW_HPP

#include "utopia_Tensor.hpp"
#include "utopia_kokkos_Base.hpp"

#include <Kokkos_View.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_Core.hpp>

namespace utopia {

    template<class KokkosView2D_>
    class MatrixView final : public Tensor<MatrixView<KokkosView2D_>, 2> {
    public:
        using KokkosView = KokkosView2D_;
        using Scalar = typename KokkosView::value_type;
        using SizeType = std::size_t;
        
        //FIXME
        typedef Kokkos::TeamPolicy<>               TeamPolicy;
        typedef Kokkos::TeamPolicy<>::member_type  MemberType;

        using Super = utopia::Tensor<MatrixView, 2>;
        using Super::Super;

        inline std::string get_class() const override
        {
            return "MatrixView";
        }

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

        UTOPIA_FUNCTION void assign(const MatrixView &other) override
        {
            UTOPIA_DEVICE_ASSERT(rows() == other.rows());
            UTOPIA_DEVICE_ASSERT(cols() == other.cols());

            const SizeType r = rows();
            const SizeType c = cols();

            Kokkos::parallel_for(
                "MatrixView::assign",
                TeamPolicy(r, Kokkos::AUTO),
                KOKKOS_LAMBDA(const MemberType &team_member) {
                    const int i = team_member.league_rank();
                    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, c), [&] (const int j) {
                    view_(i, j) = other.view_(i, j);
                });
            });
        }

        UTOPIA_FUNCTION void assign(MatrixView &&other) override
        {
            view_ = std::move(other.view_);
        }

        UTOPIA_INLINE_FUNCTION KokkosView &raw_type() { return view_; }
        UTOPIA_INLINE_FUNCTION const KokkosView &raw_type() const { return view_; }

        UTOPIA_INLINE_FUNCTION void resize(const SizeType &rows, const SizeType &cols)
        {
            Kokkos::resize(view_, rows, cols);
        }

        UTOPIA_INLINE_FUNCTION SizeType rows() const
        {
            return view_.extent(0);
        }

        UTOPIA_INLINE_FUNCTION SizeType cols() const
        {
            return view_.extent(1);
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
            KokkosBlas::scal(view_,alpha,view_);
        }

        UTOPIA_INLINE_FUNCTION void axpy(const Scalar &alpha, const MatrixView &x)
        {
            KokkosBlas::axpy(alpha,x.view_,view_);
        }

        UTOPIA_INLINE_FUNCTION Scalar dot(const MatrixView &other) const
        {
            return KokkosBlas::dot(view_, other.view_);
        }

        template<class OtherView>
        UTOPIA_INLINE_FUNCTION void multiply(const VectorView<OtherView> &right, VectorView<OtherView> &result) const
        {
            UTOPIA_DEVICE_ASSERT(!result.is_alias(right));
            UTOPIA_DEVICE_ASSERT(result.size() == rows());
            UTOPIA_DEVICE_ASSERT(right.size() == cols());

            // result.set(0.0);

            // const SizeType r = rows();
            // const SizeType c = right.size();

            // for(SizeType i = 0; i < r; ++i) {
            //     for(SizeType j = 0; j < c; ++j) {
            //         result.add(i, get(i, j) * right.get(j));
            //     }
            // }

            KokkosBlas::gemv("N", 1.0, view_, right.raw_type(), 0.0, result.raw_type());
        }


        UTOPIA_INLINE_FUNCTION void multiply(const MatrixView &right, MatrixView &result) const
        {
            UTOPIA_DEVICE_ASSERT(!result.is_alias(right));
            UTOPIA_DEVICE_ASSERT(!is_alias(result));
            UTOPIA_DEVICE_ASSERT(rows() == result.rows());
            UTOPIA_DEVICE_ASSERT(right.cols() == result.cols());
            UTOPIA_DEVICE_ASSERT(cols() == right.rows());

            const SizeType r = rows();
            const SizeType m = cols();
            const SizeType c = right.cols();

            // result.set(0.0);

            // for(SizeType i = 0; i < r; ++i) {
            //     for(SizeType j = 0; j < c; ++j) {
            //         for(SizeType k = 0; k < m; ++k) {
            //             result.add(i, j, get(i, k) * right.get(k, j));
            //         }
            //     }
            // }

            KokkosBlas::gemm("N", "N", 1.0, view_, right.view_, 0.0, result.view_);
        }

        UTOPIA_INLINE_FUNCTION void set(const Scalar &alpha)
        {
            KokkosBlas::fill(view_, alpha);
        }

        UTOPIA_FUNCTION MatrixView(const KokkosView &view)
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

        UTOPIA_INLINE_FUNCTION void wrap(const KokkosView &view)
        {
            view_ = view;
        }

        UTOPIA_INLINE_FUNCTION bool is_alias(const MatrixView &other) const
        {
            return &(view_(0, 0)) == &(other.view_(0, 0));
        }

    private:
        KokkosView view_;

        UTOPIA_FUNCTION MatrixView(const MatrixView &) {
            UTOPIA_DEVICE_ASSERT(false);
        }
    };

}

#endif //UTOPIA_KOKKOS_MATRIX_VIEW_HPP

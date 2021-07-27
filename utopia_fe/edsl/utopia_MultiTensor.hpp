#ifndef UTOPIA_FE_TENSOR_HPP
#define UTOPIA_FE_TENSOR_HPP

#include "moonolith_expanding_array.hpp"
#include "utopia_LocalFEForwardDeclarations.hpp"
#include "utopia_Traits.hpp"
#include "utopia_blas.hpp"

#include <cmath>
#include <iostream>

namespace utopia {

    static const int FEBLAS = BLAS + 1;
    template <typename T>
    class BlasMultiTensorTraits {
    public:
        using BlasTraitsT = utopia::Traits<BlasVector<T>>;
        using Scalar = T;
        using MultiScalar = utopia::MultiTensor<T, 0>;
        using Matrix = utopia::MultiTensor<utopia::BlasMatrix<T>, 2>;
        using Vector = utopia::MultiTensor<utopia::BlasVector<T>, 1>;

        using SizeType = typename BlasTraitsT::SizeType;
        using IndexSet = typename BlasTraitsT::IndexSet;
        using ScalarArray = typename BlasTraitsT::ScalarArray;
        using IndexArray = typename BlasTraitsT::IndexArray;

        enum { Backend = FEBLAS };

        static BackendInfo &backend_info() {
            static BackendInfo instance_("fe_blas");
            return instance_;
        }
    };

    using BlasMultiTensorTraitsd = utopia::BlasMultiTensorTraits<double>;
    using MultiVectord = utopia::MultiTensor<BlasVectord, 1>;
    using MultiMatrixd = utopia::MultiTensor<BlasMatrixd, 2>;
    using MultiScalard = utopia::MultiTensor<double, 0>;

    UTOPIA_MAKE_TRAITS(MultiScalard, BlasMultiTensorTraitsd, 0);
    UTOPIA_MAKE_TRAITS(MultiVectord, BlasMultiTensorTraitsd, 1);
    UTOPIA_MAKE_TRAITS_DENSE(MultiMatrixd, BlasMultiTensorTraitsd, 2);

    template <class T, int Order, int Sparsity>
    class TensorQuery<Traits<MultiTensor<T, Order>>, 0, Sparsity> {
    public:
        using Type = utopia::MultiTensor<T, 0>;
    };

    template <typename T>
    class Math<MultiTensor<T, 0>> {
    public:
        inline static MultiTensor<T, 0> abs(const MultiTensor<T, 0> &x) {
            auto n = x.size();
            MultiTensor<T, 0> ret = x;
            for (SizeType i = 0; i < n; ++i) {
                ret[i] = std::abs(ret[i]);
            }

            return ret;
        }

        inline static MultiTensor<T, 0> abs(MultiTensor<T, 0> &&x) {
            auto n = x.size();
            MultiTensor<T, 0> ret = x;
            for (SizeType i = 0; i < n; ++i) {
                ret[i] = std::abs(ret[i]);
            }

            return ret;
        }
    };

    template <typename T>
    MultiTensor<T, 0> operator*(const MultiTensor<T, 0> &left, const MultiTensor<T, 0> &right) {
        auto n = left.size();
        MultiTensor<T, 0> ret = left;
        for (SizeType i = 0; i < n; ++i) {
            ret[i] *= right[i];
        }

        return ret;
    }

    template <class Derived, int Order>
    class MultiTensorSpecificImpl {};

    template <class Derived>
    class MultiTensorSpecificImpl<Derived, 0> {
    public:
        static const int Order = 0;
        using Traits = utopia::Traits<Derived>;

        //     using SizeType     = typename Traits::SizeType;
        using Scalar = typename Traits::Scalar;

        //     // void axpy(
        //     //     const Scalar &alpha,
        //     //     const Scalar &x)
        //     // {
        //     //     auto &d = derived_();
        //     //     const Scalar offset = alpha * x;

        //     //     const SizeType n = d.size();
        //     //     for(SizeType i = 0; i < n; ++i) {
        //     //         d[i] += offset;
        //     //     }
        //     // }

        void set(const Scalar &val) {
            auto &d = derived_();
            const SizeType n = d.size();
            for (SizeType i = 0; i < n; ++i) {
                d.at(i) = val;
            }
        }

    private:
        inline Derived &derived_() { return static_cast<Derived &>(*this); }

        inline const Derived &derived_() const { return static_cast<const Derived &>(*this); }
    };

    template <class Derived>
    class MultiTensorSpecificImpl<Derived, 2> {
    public:
        static const int Order = 2;
        using Traits = utopia::Traits<Derived>;

        using T = typename Traits::Matrix;
        using SizeType = typename Traits::SizeType;
        using Scalar = typename Traits::Scalar;
        // using MultiScalar  = typename Traits::MultiScalar;

        virtual ~MultiTensorSpecificImpl() {}

        void shift_diag(const Scalar &alpha) {
            auto &d = derived_();
            const auto n = d.size();
            for (SizeType i = 0; i < n; ++i) {
                d.at(i).shift_diag(alpha);
            }
        }

        void transpose(Derived &result) const {
            auto &d = derived_();
            const auto n = d.size();

            if (result.size() != n) {
                result.resize(n);
            }

            for (SizeType i = 0; i < n; ++i) {
                d.at(i).transpose(result.at(i));
            }
        }

    private:
        inline Derived &derived_() { return static_cast<Derived &>(*this); }

        inline const Derived &derived_() const { return static_cast<const Derived &>(*this); }
    };

    template <class T_, int Order_ = Traits<T_>::Order>
    class MultiTensor final : public Tensor<MultiTensor<T_>, Order_>,
                              public MultiTensorSpecificImpl<MultiTensor<T_>, Order_> {
    public:
        using T = T_;
        static const int Order = Order_;

        using Traits = utopia::Traits<MultiTensor<T, Order>>;

        using SizeType = typename Traits::SizeType;
        using Scalar = typename Traits::Scalar;
        using MultiScalar = utopia::MultiTensor<Scalar, 0>;

        using Super = utopia::Tensor<MultiTensor, Order_>;
        using Super::Super;

        inline std::string get_class() const override { return "Multi" + Order2String<Order>::Value(); }

        MultiTensor() {}

        MultiTensor(MultiTensor &&other) : values_(std::move(other.values_)) {}

        MultiTensor(const MultiTensor &other) : values_(other.values_) {}

        MultiTensor &operator=(MultiTensor &&other) {
            values_ = std::move(other.values_);
            return *this;
        }

        MultiTensor &operator=(const MultiTensor &other) {
            if (is_alias(other)) return *this;
            values_ = other.values_;
            return *this;
        }

        template <class Expr>
        MultiTensor(const Expression<Expr> &expr) {
            // THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::construct_eval(expr.derived());
        }

        template <class Expr>
        inline MultiTensor &operator=(const Expression<Expr> &expr) {
            // THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::assign_eval(expr.derived());
            return *this;
        }

        void assign(const MultiTensor &other) override { values_ = other.values_; }

        void assign(MultiTensor &&other) override { values_ = std::move(other.values_); }

        template <class TAlpha, int OrderAlpha>
        void axpy(const MultiTensor<TAlpha, OrderAlpha> &alpha, const MultiTensor &x) {
            const auto n = size();
            assert(n == alpha.size());
            assert(n == x.size());

            for (SizeType i = 0; i < n; ++i) {
                at(i) += alpha[i] * x[i];
            }
        }

        void axpy(const Scalar &alpha, const MultiTensor &x) {
            const auto n = size();
            assert(n == x.size());

            for (SizeType i = 0; i < n; ++i) {
                at(i) += alpha * x[i];
            }
        }

        void axpy(const Scalar &alpha, const Scalar &x) {
            const Scalar offset = alpha * x;
            const SizeType n = size();
            for (SizeType i = 0; i < n; ++i) {
                at(i) += offset;
            }
        }

        void dot(const MultiTensor &other, MultiScalar &result) const {
            const auto n = size();
            assert(n == other.size());

            result.resize(n);

            for (SizeType i = 0; i < n; ++i) {
                result[i] = utopia::dot(at(i), other[i]);
            }
        }

        MultiScalar dot(const MultiTensor &other) const {
            MultiScalar result;
            dot(other, result);
            return result;
        }

        void scale(const MultiScalar &alpha) {
            const auto n = size();
            assert(n == alpha.size());

            for (SizeType i = 0; i < n; ++i) {
                at(i) *= alpha[i];
            }
        }

        template <class OtherT, int OtherOrder>
        void multiply(const MultiTensor<OtherT, OtherOrder> &other, MultiTensor<OtherT, OtherOrder> &result) const {
            const auto n = size();
            assert(n == other.size());

            if (n != result.size()) {
                result.resize(n);
            }

            for (SizeType i = 0; i < n; ++i) {
                result[i] = at(i) * other.at(i);
            }
        }

        void scale(const Scalar &alpha) {
            const auto n = size();
            for (SizeType i = 0; i < n; ++i) {
                at(i) *= alpha;
            }
        }

        template <class Op>
        void transform(const Op &op) {
            const auto n = size();
            for (SizeType i = 0; i < n; ++i) {
                transform_aux(op, at(i));
            }
        }

        template <class Op, class Derived, int TOrder>
        static void transform_aux(const Op &op, Tensor<Derived, TOrder> &t) {
            t.derived().transform(op);
        }

        template <class Op>
        static void transform_aux(const Op &op, Scalar &value) {
            value = op.apply(value);
        }

        inline SizeType size() const { return values_.size(); }

        inline void resize(const SizeType n) { values_.resize(n); }

        inline const T &operator[](const SizeType i) const {
            assert(i < size());
            return values_[i];
        }

        inline T &operator[](const SizeType i) {
            assert(i < size());
            return values_[i];
        }

        inline T &at(const SizeType i) {
            assert(i < size());
            return values_[i];
        }

        inline const T &at(const SizeType i) const {
            assert(i < size());
            return values_[i];
        }

        inline bool is_alias(const MultiTensor &other) const {
            return this == &other;
            // if we move to views this is actually better
            // return &values_[0] == &other.values_[0];
        }

        inline static bool is_alias(const Number<Scalar> &other) { return false; }

        inline void describe() const {
            const auto n = size();

            std::cout << "----------------------\n";
            for (SizeType i = 0; i < n; ++i) {
                std::cout << i << ")\n";
                disp(at(i));
            }
            std::cout << "----------------------\n";
        }

    private:
        moonolith::Storage<T> values_;
    };
}  // namespace utopia

#endif  // UTOPIA_FE_TENSOR_HPP

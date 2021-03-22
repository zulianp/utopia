#ifndef UTOPIA_LITERAL_HPP
#define UTOPIA_LITERAL_HPP

#include <cmath>
#include <iostream>
#include "utopia_BLAS_Operands.hpp"
#include "utopia_Base.hpp"
#include "utopia_Communicator.hpp"
#include "utopia_Expression.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Normed.hpp"

namespace utopia {

    template <typename T>
    class Zero {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr T value() { return static_cast<T>(0); }
    };

    template <typename T>
    class One {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr T value() { return static_cast<T>(1); }
    };

    template <typename T>
    class Math {
    public:
        inline static T abs(const T &x) { return std::abs(x); }
    };

    template <typename _Scalar>
    class Number : public Expression<Number<_Scalar>>, public BLAS1Tensor<Number<_Scalar>>, public Normed<_Scalar> {
    public:
        static const int Order = 0;
        // static const int StoreAs = UTOPIA_BY_VALUE;

        enum { StoreAs = UTOPIA_BY_VALUE };

        using Scalar = _Scalar;

        inline operator Scalar &() { return value_; }

        inline constexpr operator const Scalar &() const { return value_; }

        inline Number &operator+=(const Scalar &other) {
            value_ += other;
            return *this;
        }

        inline Number &operator-=(const Scalar &other) {
            value_ -= other;
            return *this;
        }

        inline Number &operator*=(const Scalar &other) {
            value_ *= other;
            return *this;
        }

        inline Number &operator/=(const Scalar &other) {
            value_ /= other;
            return *this;
        }

        inline bool is_alias(const Number &other) const { return this == &other; }

        inline Scalar get() const { return value_; }

        inline void set(const Scalar &value) { value_ = value; }
        inline constexpr Number(const Scalar &value = Zero<Scalar>::value()) : value_(value) {}

        inline std::string get_class() const override { return "Number"; }

        template <typename OtherScalar>
        constexpr Number(const Number<OtherScalar> &other) : value_(other.value_) {}

        template <typename OtherScalar>
        constexpr Number(Number<OtherScalar> &&other) : value_(std::move(other.value_)) {}

        template <class Derived>
        Number(const Expression<Derived> &expr) : value_(scalar_cast<Scalar>(expr.derived())) {}

        template <class Derived>
        inline void construct(const Expression<Derived> &expr) {
            value_ = scalar_cast<Scalar>(expr.derived());
        }

        template <class OtherScalar>
        inline void construct(const Number<OtherScalar> &expr) {
            value_ = scalar_cast<Scalar>(expr.value_);
        }

        template <class Any>
        inline void zeros(Any &&) {
            value_ = 0.0;
        }

        inline void construct(const Number<Scalar> &expr) { value_ = expr.value_; }

        inline void assign(const Number<Scalar> &expr) { value_ = expr.value_; }

        ///< Scalar>SWAP - swap x and y
        inline void swap(Number &x) override { std::swap(x.value_, value_); }

        ///< Scalar>SCAL - x = a*x
        inline void scale(const Scalar &a) override { value_ *= a; }

        ///< Scalar>COPY - copy x into y (this)
        inline void copy(const Number &x) override { value_ = x.value_; }

        ///< Scalar>AXPY - y = a*x + y
        inline void axpy(const Scalar &a, const Number &x) override { value_ += a * x.value_; }

        ///< Scalar>DOT - dot product
        inline Scalar dot(const Number &other) const override { return value_ * other.value_; }

        ///< Scalar>NRM2 - Euclidean norm
        inline Scalar norm2() const override { return Math<Scalar>::abs(value_); }

        ///< Scalar>ASUM - sum of absolute values
        inline Scalar norm1() const override { return Math<Scalar>::abs(value_); }

        ///< Scalar>ASUM - sum of absolute values
        inline Scalar norm_infty() const override { return Math<Scalar>::abs(value_); }

        /// I<Scalar>AMAX - index of max abs value
        inline unsigned short amax() const  // override
        {
            return 0;
        }

    private:
        Scalar value_;
    };

    template <typename T>
    Size size(const Number<T> &) {
        return {1};
    }

    template <typename T>
    inline Layout<SelfCommunicator, 1, int> layout(const Number<T> &) {
        return Layout<SelfCommunicator, 1, int>(SelfCommunicator(), 1, 1);
    }

    template <typename Scalar_>
    class Traits<Number<Scalar_>> : public Traits<Scalar_> {
    public:
        using SizeType = unsigned short;
    };

    template <typename T>
    class IsScalar {
    public:
        enum { value = false };
    };

    template <typename T>
    class IsScalar<Number<T>> {
    public:
        enum { value = true };
    };

    template <typename T>
    inline constexpr bool is_scalar() {
        return IsScalar<T>::value;
    }

    template <typename _Scalar, _Scalar Value = 0>
    class Literal : public Expression<Literal<_Scalar>> {
    public:
        using Scalar = _Scalar;

        static const int Order = 0;

        inline constexpr operator Scalar() const { return Value; }

        inline std::string get_class() const { return "Literal(TODO)"; }
    };

    template <typename T>
    void disp(const Number<T> &num, std::ostream &os = std::cout) {
        os << static_cast<T>(num) << std::endl;
    }
}  // namespace utopia
#endif  // UTOPIA_LITERAL_HPP

//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef SIMMOD_utopia_OPERATORS_HPP
#define SIMMOD_utopia_OPERATORS_HPP

#include <string>
#include <functional>
#include <cmath>

namespace utopia {

    template<class Op>
    class is_commutative {
    public:
        static const int value = 0;
    };

    template<class Op>
    class is_associative {
    public:
        static const int value = 0;
    };

    /**
     * @brief compile time information about the combination of operators
     *  true if (A Op1 B Op2 C) == (A Op1 B) Op2 == A Op1 (B Op2 C)
     */
    template<class Op1, class Op2>
    class eval_order_changable {
    public:
        static const int value = 0;
    };

    class Minus {
    public:
        std::string get_class() const { return "Minus"; }

        template<typename T>
        inline static T apply(const T &l, const T &r)
        {
            return l - r;
        }

        template<typename T>
        inline static T apply(const T &expr)
        {
            return -expr;
        }
    };

    class Plus {
    public:
        std::string get_class() const { return "Plus"; }

        template<typename T>
        inline static T apply(const T &l, const T &r)
        {
            return l + r;
        }
    };

    template<>
    class is_commutative<Plus> {
    public:
        static const int value = 1;
    };

    template<>
    class is_associative<Plus> {
    public:
        static const int value = 1;
    };

    class PlusEqual {
    public:
        std::string get_class() const { return "PlusEqual"; }

        template<typename T>
        inline static T & apply(T &l, const T &r)
        {
            return l += r;
        }
    };

    class AbsPlus {
    public:
        std::string get_class() const { return "AbsPlus"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            using std::abs;
            return abs(left) + abs(right);
        }
    };

    template<>
    class is_commutative<AbsPlus> {
    public:
        static const int value = 1;
    };

    template<>
    class is_associative<AbsPlus> {
    public:
        static const int value = 1;
    };


    class And {
    public:
        std::string get_class() const { return "And"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            return left && right;
        }
    };

    template<>
    class is_commutative<And> {
    public:
        static const int value = 1;
    };

    template<>
    class is_associative<And> {
    public:
        static const int value = 1;
    };


    class Multiplies {
    public:
        std::string get_class() const { return "Multiplies"; }

        template<typename Left, typename Right>
        inline static auto apply(const Left &left, const Right &right) -> decltype(left * right) {
            return left * right;
        }
    };

    template<>
    class is_commutative<Multiplies> {
    public:
        static const int value = 0;
    };

    template<>
    class is_associative<Multiplies> {
    public:
        static const int value = 1;
    };


    class Divides {
    public:
        std::string get_class() const { return "Divides"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            return left / right;
        }
    };

    class EMultiplies {
    public:
        std::string get_class() const { return "EMultiplies"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            return left * right;
        }
    };

    template<>
    class is_commutative<EMultiplies> {
    public:
        static const int value = 1;
    };

    template<>
    class is_associative<EMultiplies> {
    public:
        static const int value = 1;
    };


    class KroneckerProduct {
    public:
        std::string get_class() const { return "KroneckerProduct"; }
    };

    class TraceOp {
    public:
        std::string get_class() const { return "TraceOp"; }
    };

    class ApproxEqual {
    public:
        std::string get_class() const { return "ApproxEqual"; }

        template<typename T>
        inline bool operator()(const T &left, const T &right) const {
            using std::abs;
            return abs(left - right) < _tol;
        }

        ApproxEqual(const double tol = 1e-8)
                : _tol(tol) { }

        inline double tol() const {
            return _tol;
        }

        template<typename T>
        inline bool apply(const T &left, const T &right) const {
            using std::abs;
            return abs(left - right) < tol();
        }

    private:
        double _tol;
    };

    template<>
    class is_commutative<ApproxEqual> {
    public:
        static const int value = 1;
    };

    template<>
    class is_associative<ApproxEqual> {
    public:
        static const int value = 0;
    };


    //UNARY
    class Sqrt {
    public:
        std::string get_class() const { return "Sqrt"; }

        template<typename T>
        inline static T apply(const T &x) {
            using std::sqrt;
            return sqrt(x);
        }
    };

    class Pow2 {
    public:
        std::string get_class() const { return "Pow2"; }

        template<typename T>
        inline static T apply(const T &x) {
            return x * x;
        }
    };


    class Pow {
    public:
        std::string get_class() const { return "Pow"; }

        template<typename T>
        inline T apply(const T &x) const {
            return std::pow(x,a_);
        }

       Pow(const double &a):a_(a){}
       double a_;

    };

    class Log {
    public:
        std::string get_class() const { return "Log"; }

        template<typename T>
        inline static T apply(const T &x) {
            using std::log;
            return log(x);
        }
    };

    class Exp {
    public:
        std::string get_class() const { return "Exp"; }

        template<typename T>
        inline static T apply(const T &x) {
            using std::exp;
            return exp(x);
        }
    };

    class Cos {
    public:
        std::string get_class() const { return "Cos"; }

        template<typename T>
        inline static T apply(const T &x) {
            using std::cos;
            return cos(x);
        }
    };


    class Sin {
    public:
        std::string get_class() const { return "Sin"; }

        template<typename T>
        inline static T apply(const T &x) {
            using std::sin;
            return sin(x);
        }
    };

    class Transpose {
    public:
        std::string get_class() const { return "Transpose"; }
    };

    class Abs {
    public:
        std::string get_class() const { return "Abs"; }

        template<typename T>
        inline static T apply(const T &x) {
            using std::abs;
            return abs(x);
        }
    };


    template<class Expr>
    std::string GetClass() {
        return Expr().get_class();
    }

    //Special (FIND A NAME)
    template<typename  T>
    class Reciprocal {
    public:
        typedef T Scalar;

        std::string get_class() const { return "Reciprocal"; }

        template<typename T2>
        inline T2 apply(const T2 &x) const {
            return numerator_/x;
        }

        Reciprocal(const T &numerator)
                : numerator_(numerator)
        {}

        inline const T &numerator() const { return numerator_; }
    private:
        const T numerator_;
    };


    class Min {
    public:
        std::string get_class() const { return "Min"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            using std::min;
            return min(left, right);
        }
    };

    template<>
    class is_commutative<Min> {
    public:
        static const int value = 1;
    };


    class Max {
    public:
        std::string get_class() const { return "Max"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            using std::max;
            return max(left, right);
        }
    };

    template<>
    class is_commutative<Max> {
    public:
        static const int value = 1;
    };


    class IsNaNOrInf {
    public:
        std::string get_class() const { return "IsNaNOrInf"; }
        template<typename T>
        inline static T apply(const T &value) {
            return std::isnan(value) || std::isinf(value);
        }
    };

    template<typename T>
    class IsNonZero {
    public:
        std::string get_class() const { return "IsNonZero"; }

        IsNonZero(const T tol = 0.)
        : tol_(tol)
        {}

        template<typename TOther>
        inline TOther apply(const TOther &x) const {
            using std::abs;
            return abs(x) > tol_;
        }

    private:
        T tol_;
    };

    template<typename T>
    class PlusIsNonZero {
    public:
        PlusIsNonZero(const T tol = 0.)
        : is_non_zero_(tol)
        {}

        template<typename TOther>
        inline TOther apply(const TOther &x, const TOther &y) const {
            return is_non_zero_.apply(x) + is_non_zero_.apply(y);
        }

        inline const IsNonZero<T> &is_non_zero() const
        {
            return is_non_zero_;
        }

        std::string get_class() const { return "PlusIsNonZero"; }

    private:
        IsNonZero<T> is_non_zero_;

    };

    template<class T>
    class eval_order_changable<T, T> {
    public:
        static const int value = is_associative<T>::value;
    };


    template<>
    class eval_order_changable<Plus, Plus> {
    public:
        static const int value = 1;
    };

    template<>
    class eval_order_changable<Minus, Minus> {
    public:
        static const int value = 1;
    };

    template<>
    class eval_order_changable<Plus, Minus> {
    public:
        static const int value = 1;
    };

    template<>
    class eval_order_changable<Minus, Plus> {
    public:
        static const int value = 1;
    };
    
}

#endif //SIMMOD_utopia_OPERATORS_HPP

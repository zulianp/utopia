//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef SIMMOD_utopia_OPERATORS_HPP
#define SIMMOD_utopia_OPERATORS_HPP

#include <string>
#include <functional>
#include <cmath>

namespace utopia {
    class Minus {
    public:
        std::string getClass() const { return "Minus"; }

        template<typename T>
        inline static T apply(const T &l, const T &r)
        {
            return l - r;
        }
    };

    class Plus {
    public:
        std::string getClass() const { return "Plus"; }

        template<typename T>
        inline static T apply(const T &l, const T &r)
        {
            return l + r;
        }
    };

    class PlusEqual {
    public:
        std::string getClass() const { return "PlusEqual"; }

        template<typename T>
        inline static T & apply(T &l, const T &r)
        {
            return l += r;
        }
    };

    class AbsPlus {
    public:
        std::string getClass() const { return "AbsPlus"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            using std::abs;
            return abs(left) + abs(right);
        }
    };

    class And {
    public:
        std::string getClass() const { return "And"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            return left && right;
        }
    };

    class Multiplies {
    public:
        std::string getClass() const { return "Multiplies"; }

        template<typename Left, typename Right>
        inline static auto apply(const Left &left, const Right &right) -> decltype(left * right) {
            return left * right;
        }
    };

    class Divides {
    public:
        std::string getClass() const { return "Divides"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            return left / right;
        }
    };

    class EMultiplies {
    public:
        std::string getClass() const { return "EMultiplies"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            return left * right;
        }
    };

    class KroneckerProduct {
    public:
        std::string getClass() const { return "KroneckerProduct"; }
    };

    class TraceOp {
    public:
        std::string getClass() const { return "TraceOp"; }
    };

    class ApproxEqual {
    public:
        std::string getClass() const { return "ApproxEqual"; }

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

    //UNARY
    class Sqrt {
    public:
        std::string getClass() const { return "Sqrt"; }

        template<typename T>
        inline static T apply(const T &x) {
            using std::sqrt;
            return sqrt(x);
        }
    };

    class Pow2 {
    public:
        std::string getClass() const { return "Pow2"; }

        template<typename T>
        inline static T apply(const T &x) {
            return x * x;
        }
    };

    class Transpose {
    public:
        std::string getClass() const { return "Transpose"; }
    };

    class Abs {
    public:
        std::string getClass() const { return "Abs"; }

        template<typename T>
        inline static T apply(const T &x) {
            using std::abs;
            return abs(x);
        }
    };


    template<class Expr>
    std::string GetClass() {
        return Expr().getClass();
    }

    //Special (FIND A NAME)
    template<typename  T>
    class Reciprocal {
    public:
        typedef T Scalar;
        
        std::string getClass() const { return "Reciprocal"; }

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
        std::string getClass() const { return "Min"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            using std::min;
            
            return min(left, right);
        }
    };

    class Max {
    public:
        std::string getClass() const { return "Max"; }

        template<typename T>
        inline static T apply(const T &left, const T &right) {
            using std::max;
            return max(left, right);
        }
    };
}

#endif //SIMMOD_utopia_OPERATORS_HPP

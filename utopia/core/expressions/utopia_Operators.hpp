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
        inline static std::minus<T> Fun() {
            return std::minus<T>();
        }
    };

    class Plus {
    public:
        std::string getClass() const { return "Plus"; }

        template<typename T>
        inline static std::plus<T> Fun() {
            return std::plus<T>();
        }
    };

    class AbsPlus {
    public:
        std::string getClass() const { return "AbsPlus"; }

        template<typename T>
        class Function {
        public:
            inline T operator()(const T &left, const T &right) const {
                using std::abs;
                return abs(left) + abs(right);
            }
        };

        template<typename T>
        inline static Function<T> Fun() {
            return Function<T>();
        }

    };

    class And {
    public:
        std::string getClass() const { return "And"; }

        template<typename T>
        class Function {
        public:
            inline T operator()(const T &left, const T &right) const {
                return left && right;
            }
        };

        template<typename T>
        inline static Function<T> Fun() {
            return Function<T>();
        }

    };


    class Multiplies {
    public:
        std::string getClass() const { return "Multiplies"; }

        template<typename T>
        inline static std::multiplies<T> Fun() {
            return std::multiplies<T>();
        }
    };

    class Divides {
    public:
        std::string getClass() const { return "Divides"; }

        template<typename T>
        inline static std::divides<T> Fun() {
            return std::divides<T>();
        }
    };



    class EMultiplies {
    public:
        std::string getClass() const { return "EMultiplies"; }

        template<typename T>
        inline static std::multiplies<T> Fun() {
            return std::multiplies<T>();
        }
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

        double getTol() const {
            return _tol;
        }

    private:
        double _tol;
    };

    //UNARY
    class Sqrt {
    public:
        std::string getClass() const { return "Sqrt"; }

        template<typename T>
        class Function {
        public:
            inline T operator()(const T &x) const {
                using std::sqrt;
                return sqrt(x);
            }
        };

        template<typename T>
        inline static Function<T> Fun() {
            return Function<T>();
        }
    };

    class Pow2 {
    public:
        std::string getClass() const { return "Pow2"; }

        template<typename T>
        class Function {
        public:
            inline T operator()(const T &x) const {
                return x * x;
            }
        };

        template<typename T>
        inline static Function<T> Fun() {
            return Function<T>();
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
        class Function {
        public:
            inline T operator()(const T &x) const {
                using std::abs;
                return abs(x);
            }
        };

        template<typename T>
        inline static Function<T> Fun() {
            return Function<T>();
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
        class Function {
        public:
            inline T operator()(const T &left, const T &right) const {
                using std::min;
                return min(left, right);
            }
        };

        template<typename T>
        inline static Function<T> Fun() {
            return Function<T>();
        }
    };

    class Max {
    public:
        std::string getClass() const { return "Max"; }

        template<typename T>
        class Function {
        public:
            inline T operator()(const T &left, const T &right) const {
                using std::max;
                return max(left, right);
            }
        };

        template<typename T>
        inline static Function<T> Fun() {
            return Function<T>();
        }
    };
}

#endif //SIMMOD_utopia_OPERATORS_HPP

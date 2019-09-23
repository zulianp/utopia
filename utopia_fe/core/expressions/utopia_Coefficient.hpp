#ifndef UTOPIA_FE_COEFFICIENT_HPP
#define UTOPIA_FE_COEFFICIENT_HPP

#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include <utility>
#include <functional>

namespace utopia {

    template<typename T, int Order_>
    class ConstantCoefficient final : public Expression< ConstantCoefficient<T, Order_> > {
    public:
        static const int Order = Order_;

        typedef typename Traits<T>::Scalar Scalar;

        ConstantCoefficient(const T &value)
        : value_(value)
        {}


        inline const T &expr() const
        {
            return value_;
        }

        inline operator const T&() const
        {
            return value_;
        }

        inline std::string get_class() const override
        {
            return "ConstantCoefficient" + std::to_string(Order);
        }


        //always returns value_
        inline const T &operator[](const int) const
        {
            return value_;
        }

        template<class Point, class Result>
        void eval(const Point &, Result &result) const {
            result = value_;
        }

    private:
        T value_;
    };

    template<typename T>
    class ConstantCoefficient<T, 0> final : public Expression< ConstantCoefficient<T, 0> > {
    public:

        enum {
            Order = 0
        };

        // typedef T Scalar;
        using Scalar = typename Traits<T>::Scalar;

        ConstantCoefficient(const T &value)
        : value_(value)
        {}


        inline const T &expr() const
        {
            return value_;
        }

        inline operator const T&() const
        {
            return value_;
        }

        //always returns value_
        inline const T &operator[](const int) const
        {
            return value_;
        }

        template<class Point, class Result>
        void eval(const Point &, Result &result) const {
            result = value_;
        }

        inline std::string get_class() const override
        {
            return "ConstantCoefficient0";
        }
    private:
        T value_;
    };

    template<class T>
    inline ConstantCoefficient<T, 0> coeff(const T &expr) {
        return ConstantCoefficient<T, 0>(expr);
    }

    template<class T, int Order>
    inline ConstantCoefficient<Tensor<T, Order> , Order> coeff(const Tensor<T, Order> &expr) {
        return ConstantCoefficient<Tensor<T, Order>, Order>(expr);
    }


    template<class T>
    inline ConstantCoefficient<T, 1> vec_coeff(const T &expr) {
        return ConstantCoefficient<T, 1>(expr);
    }


    template<class T>
    inline ConstantCoefficient<T, 2> mat_coeff(const T &expr) {
        return ConstantCoefficient<T, 2>(expr);
    }




    template<class T>
    class Traits< ConstantCoefficient<T, 0> > {
    public:
        enum {
            FILL_TYPE = utopia::FillType::DENSE
        };

        enum {
            Order = 0
        };

        typedef T Scalar;
    };


    template<class T>
    class Traits< ConstantCoefficient<T, 1> > {
    public:
        enum {
            FILL_TYPE = utopia::FillType::DENSE
        };

        enum {
            Order = 1
        };

        typedef typename Traits<T>::Scalar Scalar;
    };


    template<class T>
    class Traits< ConstantCoefficient<T, 2> > {
    public:
        enum {
            FILL_TYPE = utopia::FillType::DENSE
        };

        enum {
            Order = 2
        };

        typedef typename Traits<T>::Scalar Scalar;
    };


    template<class Fun, typename T, int Order_>
    class FunctionCoefficient final : public Expression< FunctionCoefficient<Fun, T, Order_> > {
    public:
            enum {
                Order = Order_
            };

            typedef typename Traits<T>::Scalar Scalar;

            FunctionCoefficient(const Fun &fun)
            : fun_(fun)
            {}


            template<class Point, class Result>
            void eval(const Point &p, Result &result) const {
                fun_(p, result);
            }

            inline std::string get_class() const override
            {
                return "FunctionCoefficient" + std::to_string(Order);
            }

            inline const Fun &fun() const
            {
                return fun_;
            }

            inline Fun &fun()
            {
                return fun_;
            }

        private:
            Fun fun_;
    };

    template<class Fun, typename T>
    class FunctionCoefficient<Fun, T, 0> final : public Expression< FunctionCoefficient<Fun, T, 0> > {
    public:
            enum {
                Order = 0
            };

            typedef T Scalar;

            FunctionCoefficient() = default;

            FunctionCoefficient(const Fun &fun)
            : fun_(fun)
            {}


            template<class Point, class Result>
            void eval(const Point &p, Result &result) const {
                result = fun_(p);
            }


            inline std::string get_class() const override
            {
                return "FunctionCoefficient0";
            }

            inline const Fun &fun() const
            {
                return fun_;
            }

            inline Fun &fun()
            {
                return fun_;
            }

        private:
            Fun fun_;
    };

    template<class Fun, class T, int Order_>
    class Traits< FunctionCoefficient<Fun, T, Order_> > {
    public:
        enum {
            FILL_TYPE = utopia::FillType::DENSE
        };

        enum {
            Order = Order_
        };

        typedef typename Traits<T>::Scalar Scalar;
    };

    template<class Fun, class T>
    class Traits< FunctionCoefficient<Fun, T, 0> > {
    public:
        enum {
            FILL_TYPE = utopia::FillType::DENSE
        };

        enum {
            Order = 0
        };

        typedef T Scalar;
    };

    template<class Point, class T>
    inline FunctionCoefficient<T (*)(const Point&), T, 0> coeff(T (*fun)(const Point&)) {
        return fun;
    }

    template<class Point, class T>
    inline FunctionCoefficient<std::function<T(const Point&)> , T, 0> coeff(std::function<T(const Point&)> fun) {
        return fun;
    }



    template<class T, int N>
    inline FunctionCoefficient<std::function<std::array<T, N>(const std::array<T, N> &)>, T, 1> coeff(std::function<std::array<T, N>(const std::array<T, N> &)> fun) {
        return fun;
    }

    template<class Point, class T>
    inline FunctionCoefficient<void (*)(const Point&, T &), T, 1> vec_coeff(void (*fun)(const Point&, T &)) {
        return fun;
    }

    template<class Point, class T>
    inline FunctionCoefficient<std::function<void(const Point&, T &)> , T, 1> vec_coeff(std::function<void(const Point&, T&)> fun) {
        return fun;
    }

    template<class Point, class T>
    inline FunctionCoefficient<void (*)(const Point&, T &), T, 2> mat_coeff(void (*fun)(const Point&, T &)) {
        return fun;
    }
}

#endif //UTOPIA_FE_COEFFICIENT_HPP

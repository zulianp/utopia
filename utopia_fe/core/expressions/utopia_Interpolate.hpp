#ifndef UTOPIA_FE_INTERPOLATE_HPP
#define UTOPIA_FE_INTERPOLATE_HPP

#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_Any.hpp"

#include <utility>

namespace utopia {

    template<class Coefficient, class Fun>
    class Interpolate : public Expression< Interpolate<Coefficient, Fun> > {
    public:

        typedef typename Traits<Coefficient>::Scalar Scalar;

        static const int Order = Traits<Coefficient>::Order;


        template<class Coeff>
        Interpolate(Coeff &&coeff, Fun &fun)
        : coeff_(std::forward<Coeff>(coeff)), fun_(fun)
        {}

        inline Fun &fun()
        {
            return fun_;
        }

        inline const Fun &fun() const
        {
            return fun_;
        }

        inline Coefficient &coefficient()
        {
            return coeff_;
        }

        inline const Coefficient &coefficient() const
        {
            return coeff_;
        }

        inline std::string get_class() const override {
            return "Interpolate";
        }

    private:
        UTOPIA_STORE(Coefficient) coeff_;
        UTOPIA_STORE(Fun) fun_;
    };


    template<class Coefficient, class Fun>
    class Interpolate<const Coefficient, Fun> : public Expression< Interpolate<const Coefficient, Fun> > {
    public:

        typedef typename Traits<Coefficient>::Scalar Scalar;

        static const int Order = Traits<Coefficient>::Order;

        template<class Coeff>
        Interpolate(Coeff &&coeff, Fun &fun)
        : coeff_(std::forward<Coeff>(coeff)), fun_(fun)
        {}

        inline Fun &fun()
        {
            return fun_;
        }

        inline const Fun &fun() const
        {
            return fun_;
        }

        inline const Coefficient &coefficient() const
        {
            return coeff_;
        }

        inline std::string get_class() const override {
            return "Interpolate";
        }

    private:
        UTOPIA_STORE_CONST(Coefficient) coeff_;
        UTOPIA_STORE(Fun) fun_;
    };

    template<class Coefficient, class Fun>
    using GradInterpolate = Gradient<Interpolate<Coefficient, Fun> >;


    template<>
    class Interpolate<Any, Any> {};

    template<class Coefficient, class Fun>
    class Traits< Interpolate<Coefficient, Fun> > : public Traits<Fun> {
    public:
        static const int FILL_TYPE = FillType::DENSE;
    };


    template<class Coefficient, class Fun>
    inline Interpolate<Coefficient, Fun> interpolate(Coefficient &coeff, Fun &fun)
    {
        return Interpolate<Coefficient, Fun>(coeff, fun);
    }

    template<class Coefficient, class Fun>
    inline Interpolate<const Coefficient, Fun> interpolate(const Coefficient &coeff, Fun &fun)
    {
        return Interpolate<const Coefficient, Fun>(coeff, fun);
    }


    template<class Coefficient, class Fun, class Tensor>
    inline Interpolate<Coefficient, Fun> interpolate(Coefficient &coeff, Fun &fun, Tensor &&tensor)
    {
        return Interpolate<Coefficient, Fun>(coeff, fun, tensor);
    }

    //FIXME: @deprecated remove once refactor is finished
    template<class T, int Order, class Fun, class Tensor>
    inline Interpolate<ConstantCoefficient<T, Order>, Fun> interpolate(const ConstantCoefficient<T, Order> &coeff, Fun &fun, Tensor &&tensor)
    {
        return Interpolate<ConstantCoefficient<T, Order>, Fun>(coeff, fun, tensor);
    }

}

#endif //UTOPIA_FE_INTERPOLATE_HPP

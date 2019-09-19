#ifndef UTOPIA_VALUE_FUNCTION_HPP
#define UTOPIA_VALUE_FUNCTION_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_StoreAs.hpp"
#include <functional>
#include <memory>

namespace utopia {

    template<class Expr, class Out_, int Order_, class Fun>
    class Filter : public Expression< Filter<Expr, Out_, Order_, Fun> > {
    public:
        static const int Order = Order_;
        using Type = Out_;
        typedef double Scalar;

        Filter(const Expr &expr, Fun fun)
        : expr_(expr), fun_(fun)
        {}

        std::string get_class() const override { return "Filter<" + expr_.get_class() + ">"; }

        inline const Expr &expr() const
        {
            return expr_;
        }

        template<typename In>
        inline Out_ eval(const In &in) const {
            return fun_(in);
        }

        template<typename In>
        inline std::vector<Out_> eval(const std::vector<In> &in) const {
            auto n = in.size();
            std::vector<Out_> ret(n);

            for(std::size_t i = 0; i < n; ++i) {
                ret[i] = fun_(in[i]);
            }

            return ret;
        }

        template<typename In>
        inline std::vector<std::vector<Out_>> eval(const std::vector<std::vector<In>> &in) const {
            auto n = in.size();
            std::vector<std::vector<Out_>> ret(n);

            for(std::size_t i = 0; i < n; ++i) {
                auto m = in[i].size();
                ret[i].resize(m);
                for(std::size_t j = 0; j < m; ++j) {
                    ret[i][j] = fun_(in[i][j]);
                }
            }

            return ret;
        }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
        Fun fun_;
    };

    template<class Expr, class Out, int Order, class Fun>
    class Traits<Filter<Expr, Out, Order, Fun>> : public Traits<Out> {};

    template<class Expr, class In, class Out, int Order = Traits<Out>::Order>
    Filter<Expr, Out, Order, std::function<Out(In&&)>> filter(
        const Expr &expr, std::function<Out(In&&)> fun)
    {
        return Filter<Expr, Out, Order, std::function<Out(In&&)>>(expr, fun);
    }

    template<class Expr, class In, class Out, int Order = Traits<Out>::Order>
    Filter<Expr, Out, Order, std::function<Out(const In&)>> filter(
        const Expr &expr, std::function<Out(const In&)> fun)
    {
        return Filter<Expr, Out, Order, std::function<Out(const In&)>>(expr, fun);
    }

    template<class Expr, class In, class Out, int Order = Traits<Out>::Order>
    Filter<Expr, Out, Order, Out(*)(const In&)> filter(
        const Expr &expr, Out(*fun)(const In&))
    {
        return Filter<Expr, Out, Order, Out(*)(const In&)>(expr, fun);
    }
}

#endif //UTOPIA_VALUE_FUNCTION_HPP

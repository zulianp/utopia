#ifndef UTOPIA_EQUATIONS_HPP
#define UTOPIA_EQUATIONS_HPP

#include <utility>

namespace utopia {

    template<class Eqs, int Begin, int End>
    class EquationIterator {
    public:
        template<class Fun>
        void visit(Fun fun)
        {
            fun(Begin, eqs.template get<Begin>());
            EquationIterator<Eqs, Begin + 1, End> next(eqs);
            next.visit(fun);
        }

        EquationIterator(Eqs &eqs) : eqs(eqs) {}
        Eqs &eqs;
    };

    template<class Eqs, int Begin>
    class EquationIterator<Eqs, Begin, Begin> {
    public:
        template<class Fun>
        void visit(Fun) {}

        EquationIterator(const Eqs &){}
    };

    template<class... Equation>
    class Equations : public Expression< Equations<Equation...> > {
    public:

        static const int n_equations = std::tuple_size< std::tuple<Equation...> >::value;

        Equations(const Equation &...eqs)
        : eqs_(eqs...)
        { }

        template<int Index>
        inline auto get() const -> const typename std::tuple_element<Index, std::tuple<Equation...>>::type
        {
            return std::get<Index>(eqs_);
        }

        template<int Index>
        inline auto get() -> typename std::tuple_element<Index, std::tuple<Equation...>>::type
        {
            return std::get<Index>(eqs_);
        }

        template<class Fun>
        void each(Fun fun)
        {
            EquationIterator<Equations, 0, n_equations> iter(*this);
            iter.visit(fun);
        }

        template<class Fun>
        void each(Fun fun) const
        {
            EquationIterator<const Equations, 0, n_equations> iter(*this);
            iter.visit(fun);
        }


        inline std::string get_class() const override
        {
            return "Equations";
        }

    private:
        std::tuple<Equation...> eqs_;
    };


    template<class... Equation>
    inline Equations<Equation...> equations(const Equation &... eqs)
    {
        return Equations<Equation...>(eqs...);
    }

    template<class First, class...Rest>
    class GetFirst {
    public:
        typedef First Type;
    };

    template<class... Eqs, class AssemblyContext>
    class FunctionalTraits<Equations<Eqs...>, AssemblyContext>  {
    public:

        class FunctionalType {
        public:
            template<class Eq>
            void operator()(const int index, const Eq &eq)
            {
                type = std::max(type, FunctionalTraits<Eq, AssemblyContext>::type(eq, ctx));
            }

            FunctionalType(const AssemblyContext &ctx) : ctx(ctx), type(0) {}

            const AssemblyContext &ctx;
            int type;
        };

        class FunctionalOrder {
        public:

            template<class Eq>
            void operator()(const int index, const Eq &eq)
            {
                order = std::max(order,  FunctionalTraits<Eq, AssemblyContext>::order(eq, ctx));
            }

            FunctionalOrder(const AssemblyContext &ctx) : ctx(ctx), order(0) {}

            const AssemblyContext &ctx;
            int order;
        };

        inline static int type(const Equations<Eqs...> &expr,  const AssemblyContext &ctx)
        {
            FunctionalType ft(ctx);
            expr.template each<FunctionalType &>(ft);
            return ft.type;
        }

        inline static int order(const Equations<Eqs...> &expr, const AssemblyContext &ctx)
        {
            FunctionalOrder fo(ctx);
            expr.template each<FunctionalOrder &>(fo);
            return fo.order;
        }
    };
}

#endif //UTOPIA_EQUATIONS_HPP

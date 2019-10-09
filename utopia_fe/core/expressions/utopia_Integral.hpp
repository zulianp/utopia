#ifndef UTOPIA_INTEGRAL_HPP
#define UTOPIA_INTEGRAL_HPP

#include "utopia_Expression.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_FEExpression.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_FunctionalTraits.hpp"

namespace utopia {



    template<class Expr_>
    class Integral : public Expression< Integral<Expr_> >, public FEExpression {
    public:
        typedef Expr_ Expr;
        static const int Order = Expr::Order;

        enum Type {
            SURFACE = 0,
            VOLUME = 1,
            INNER_VOLUME_ONLY = 2
        };

        typedef typename Expr::Scalar Scalar;

        std::string get_class() const override { return "Integral<" + expr_.get_class() + ">"; }

        Integral(const Expr &expr, const int block_id = -1, const Type type = VOLUME)
        : expr_(expr), block_id_(block_id), integral_id_(-1), type_(type)
        {}

        inline const Expr &expr() const
        {
            return expr_;
        }

        inline int block_id() const
        {
            return block_id_;
        }

        inline int has_block_id() const
        {
            return block_id_ != -1;
        }

        inline bool has_integral_id() const
        {
            return integral_id_ != -1;
        }

        inline void set_integral_id(const int id)
        {
            integral_id_ = id;
        }

        inline bool is_surface() const
        {
            return (type_ == SURFACE);
        }

        inline bool is_volume() const
        {
            return (type_ == VOLUME);
        }

        inline bool is_inner_volume_only() const
        {
            return (type_ == INNER_VOLUME_ONLY);
        }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
        int block_id_;
        int integral_id_;
        Type type_;
    };



    template<class Derived>
    inline Integral<Derived> integral(const Expression<Derived> &expr) {
        static_assert(!IsSubTree<Integral<utopia::Any>, Derived>::value, "nested integrals are not allowed");
        return Integral<Derived>(expr.derived());
    }

    template<class Derived>
    inline Integral<Derived> integral(const Expression<Derived> &expr, const int block_id) {
        static_assert(!IsSubTree<Integral<utopia::Any>, Derived>::value, "nested integrals are not allowed");
        return Integral<Derived>(expr.derived(), block_id);
    }

    template<class Expr>
    class Traits< Integral<Expr> > : public Traits<Expr> {
    public:
        enum {
            FILL_TYPE = utopia::FillType::DENSE
        };
    };

    template<class Expr, class AssemblyContext>
    class FunctionalTraits<Integral<Expr>, AssemblyContext>  {
    public:
        inline static int type(const Integral<Expr> &expr,  const AssemblyContext &ctx) { return FunctionalTraits<Expr, AssemblyContext>::type(expr.expr(), ctx);  }
        inline static int order(const Integral<Expr> &expr, const AssemblyContext &ctx) { return FunctionalTraits<Expr, AssemblyContext>::order(expr.expr(), ctx); }
    };


    class Differential {
    public:
        constexpr Differential(const int block_id = -1) noexcept : block_id(block_id) {}
        const int block_id;
    };

    class SurfaceDifferential {
    public:
        constexpr SurfaceDifferential(const int side_set_id = -1) noexcept : side_set_id(side_set_id) {}
        const int side_set_id;
    };

    template<class Derived>
    inline Integral<Derived> surface_integral(const Expression<Derived> &expr, const int side_set_id = -1) {
        static_assert(!IsSubTree<Integral<utopia::Any>, Derived>::value, "nested integrals are not allowed");
        return Integral<Derived>(expr.derived(), side_set_id, Integral<Derived>::SURFACE);
    }

    template<class Derived>
    inline Integral<Derived> inner_volume_only_integral(const Expression<Derived> &expr, const int side_set_id = -1) {
        static_assert(!IsSubTree<Integral<utopia::Any>, Derived>::value, "nested integrals are not allowed");
        return Integral<Derived>(expr.derived(), side_set_id, Integral<Derived>::INNER_VOLUME_ONLY);
    }

     // inline constexpr Differential dV(const int block_id = -1)
     // {
     // 	return Differential(block_id);
     // }

    static const Differential dX;
    static const SurfaceDifferential dS;

     template<class Derived>
     inline Integral<Derived> operator *(const Expression<Derived> &expr, const Differential &d) {
         static_assert(!IsSubTree<Integral<utopia::Any>, Derived>::value, "nested integrals are not allowed");
         return Integral<Derived>(expr.derived(), d.block_id);
     }

     template<class Derived>
     inline Integral<Derived> operator *(const Expression<Derived> &expr, const SurfaceDifferential &d) {
         static_assert(!IsSubTree<Integral<utopia::Any>, Derived>::value, "nested integrals are not allowed");
         return Integral<Derived>(expr.derived(), d.side_set_id, Integral<Derived>::SURFACE);
     }
}

#endif //UTOPIA_INTEGRAL_HPP

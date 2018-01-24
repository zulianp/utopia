#ifndef UTOPIA_INTEGRAL_HPP
#define UTOPIA_INTEGRAL_HPP 

#include "utopia_Expression.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_FEExpression.hpp"
#include "utopia_FEIsSubTree.hpp"

namespace utopia {

	template<class Expr_>
	class Integral : public Expression< Integral<Expr_> >, public FEExpression {
	public:
		typedef Expr_ Expr;
		static const int Order = Expr::Order;


		typedef typename Expr::Scalar Scalar;

		std::string getClass() const override { return "Integral<" + expr_.getClass() + ">"; }

		Integral(const Expr &expr, const int block_id = -1, const bool is_surface = false)
		: expr_(expr), block_id_(block_id), integral_id_(-1), is_surface_(is_surface)
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
			return is_surface_;
		}

		inline bool is_volume() const
		{
			return !is_surface_;
		}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
		int block_id_;
		int integral_id_;
		bool is_surface_;
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
	inline Integral<Derived> surface_integral(const Expression<Derived> &expr, const int side_set_id) {
		static_assert(!IsSubTree<Integral<utopia::Any>, Derived>::value, "nested integrals are not allowed");
		return Integral<Derived>(expr.derived(), side_set_id, true);
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
	 	return Integral<Derived>(expr.derived(), d.side_set_id, true);
	 }
}

#endif //UTOPIA_INTEGRAL_HPP

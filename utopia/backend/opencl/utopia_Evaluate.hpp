#ifndef UTOPIA_EVALUATE_HPP
#define UTOPIA_EVALUATE_HPP 

#include "utopia_Base.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_CLTraits.hpp"

namespace utopia {
	template<class Expr, int _Order = Expr::Order>
	class Evaluate : public Expression< Evaluate<Expr, _Order> > {
	public:
		typedef void Type;
		typedef typename Expr::Scalar Scalar;
		typedef EXPR_TYPE(CLTraits<Scalar>, Expr) TensorT;

		enum {
			Order = _Order
		};

		enum {
		    StoreAs = UTOPIA_DEFAULT_EXPRESSION_STORAGE
		};

		Evaluate(const Expr &expr)
		: expr_(expr), backend_tensor_(std::make_shared<TensorT>())
		{}

		std::string getClass() const
		{
			return "Evaluate<" + expr_.getClass() + ">";
		}

		inline const Expr &expr() const
		{
			return expr_;
		}

		inline TensorT &backend_tensor()
		{
			return *backend_tensor_;
		}

		inline TensorT &backend_tensor() const
		{
			return *backend_tensor_;
		}

		inline const std::shared_ptr<TensorT> &backend_tensor_ptr() const
		{
			return backend_tensor_;
		}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
		std::shared_ptr<TensorT> backend_tensor_;
	};


	template<class Expr>
	class Evaluate<Expr, 0> {
	public:
		typedef void Type;
		typedef typename Expr::Scalar Scalar;
		typedef Scalar TensorT;

		enum {
			Order = 0
		};

		enum {
		    StoreAs = UTOPIA_DEFAULT_EXPRESSION_STORAGE
		};

		Evaluate(const Expr &expr)
		: expr_(expr), value_(std::make_shared<Scalar>(0))
		{}

		std::string getClass() const
		{
			return "Evaluate<" + expr_.getClass() + ">";
		}

		inline const Expr &expr() const
		{
			return expr_;
		}

		inline  void set_value(const Scalar value) const
		{
			*value_ = value;
		}


		inline Scalar get_value() const
		{
			return *value_;
		}


	private:
		UTOPIA_STORE_CONST(Expr) expr_;
		std::shared_ptr<Scalar> value_;
	};

	template<class Expr>
	Evaluate<Expr> make_evaluate(const Expr &expr)
	{
		return expr;
	}


	template<class Tensor, int Order>
	const Wrapper<Tensor, Order> &make_evaluate(const Wrapper<Tensor, Order> &tensor)
	{
		return tensor;
	}

	template<class Tensor, int Order>
	Wrapper<Tensor, Order> && make_evaluate(Wrapper<Tensor, Order> &&tensor)
	{
		return tensor;
	}


	template<class Expr>
	cl::Buffer &get_buffer(const Evaluate<Expr> &v)
	{
		return v.backend_tensor().buffer;
	}

	// template<class Expr>
	// cl::Buffer &get_buffer(const Evaluate<Expr, 0> &v)
	// {
	// 	return v.buffer();
	// }

	template<class Expr>
	auto size(const Evaluate<Expr> &expr) -> decltype(size(expr.expr()))
	{
		return size(expr.expr());
	}


	template<class Expr, int Order>
	class Traits< Evaluate<Expr, Order> > : public Traits<Expr> {};


	template<class InnerExpr>
	class TreeProperties< Evaluate<InnerExpr> > {
	public:
		enum { greatest_tensor_order = Evaluate<InnerExpr>::Order };
		enum { smallest_tensor_order = Evaluate<InnerExpr>::Order  };
		enum { has_mat_mat_mul 		 = TreeProperties<InnerExpr>::has_mat_mat_mul 		};
		enum { has_differentiable_sub_tree = TreeProperties<InnerExpr>::has_differentiable_sub_tree  };
	};
}

#endif //UTOPIA_EVALUATE_HPP


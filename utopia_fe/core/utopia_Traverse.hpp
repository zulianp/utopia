#ifndef UTOPIA_TRAVERSE_HPP
#define UTOPIA_TRAVERSE_HPP 

#include "utopia_FEForwardDeclarations.hpp"
#include <iostream>

namespace utopia {
	static const int TRAVERSE_CONTINUE = 0;
	static const int TRAVERSE_STOP = 1;
	static const int TRAVERSE_SKIP_SUBTREE = 2;

	template<class Expr, class Visitor>
	inline static int traverse(const Expr &expr, Visitor &visitor)
	{
		std::cout << "[Error] Encountered unhandled expression: " << expr.getClass() << std::endl;
		return TRAVERSE_CONTINUE;
	}

	template<class Expr, class Visitor>
	inline static int traverse(const Integral<Expr> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr))
		{
			case TRAVERSE_CONTINUE:
			{
				return traverse(expr.expr(), visitor);
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
	}

	template<class Expr, class Visitor>
	inline static int traverse(const Negate<Expr> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr)) {
			case TRAVERSE_CONTINUE:
			{
				return traverse(expr.expr(), visitor);
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
	}
	
	
	template<class Expr, class Operation, class Visitor>
	inline static int traverse(const Unary<Expr, Operation> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr)) {
			case TRAVERSE_CONTINUE:
			{
				return traverse(expr.expr(), visitor);
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
	}
	
	template<class Expr, class Operation, class Visitor>
	inline static int traverse(const Reduce<Expr, Operation> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr))
		{
			case TRAVERSE_CONTINUE:
			{
				return traverse(expr.expr(), visitor);
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
	}
	
	template<class Expr, class Visitor>
	inline static int traverse(const Divergence<Expr> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr))
		{
			case TRAVERSE_CONTINUE:
			{
				return traverse(expr.expr(), visitor);
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
	}

	template<class Expr, class Visitor>
	inline static int traverse(const Transposed<Expr> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr))
		{
			case TRAVERSE_CONTINUE:
			{
				return traverse(expr.expr(), visitor);
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
	}
	
	template<class Expr, class Visitor>
	inline static int traverse(const Gradient<Expr> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr))
		{
			case TRAVERSE_CONTINUE:
			{
				return traverse(expr.expr(), visitor);
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
	}
	
	template<class Expr, class Visitor>
	inline static int traverse(const Curl<Expr> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr))
		{
			case TRAVERSE_CONTINUE:
			{
				return traverse(expr.expr(), visitor);
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}	
	}
	
	template<class Left, class Right, class Operation, class Visitor>
	inline static int traverse(const Binary<Left, Right, Operation> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr))
		{
			case TRAVERSE_CONTINUE:
			{
				if(traverse(expr.left(), visitor) == TRAVERSE_CONTINUE) {
					return traverse(expr.right(), visitor);
				}
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
		
	}
	
	template<class Left, class Right, class Operation, class Visitor>
	inline static int traverse(const Binary<Left, Number<Right>, Operation> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr))
		{
			case TRAVERSE_CONTINUE:
			{
				return traverse(expr.left(), visitor);
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
	}
	
	template<class Left, class Right, class Operation, class Visitor>
	inline static int traverse(const Binary<Left, BlockVar<Right>, Operation> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr))
		{
			case TRAVERSE_CONTINUE:
			{
				return traverse(expr.left(), visitor);
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
	}
	
	template<class Left, class Right, class Visitor>
	inline static int traverse(const Multiply<Left, Right> &expr, Visitor &visitor)
	{
		switch(visitor.visit(expr))
		{
			case TRAVERSE_CONTINUE:
			{
				if(traverse(expr.left(), visitor) == TRAVERSE_CONTINUE) {
					return traverse(expr.right(), visitor);
				}
			}

			case TRAVERSE_STOP:
			{
				return TRAVERSE_STOP;
			}

			case TRAVERSE_SKIP_SUBTREE:
			{
				return TRAVERSE_CONTINUE;
			}

			default: {
				std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
				return TRAVERSE_STOP;
			}
		}
	}
	
	template<class FunctionSpaceT, class Visitor>
	inline static int traverse(const TrialFunction<FunctionSpaceT> &expr, Visitor &visitor)
	{
		return visitor.visit(expr);
	}

	template<class FunctionSpaceT, class Visitor>
	inline static int traverse(const TestFunction<FunctionSpaceT> &expr, Visitor &visitor)
	{
		return visitor.visit(expr);
	}

	template<class Expr>
	class FindExpression {
	public:

		template<class Any>
		inline constexpr static int visit(const Any &) { return TRAVERSE_CONTINUE; }
		
		int visit(const Expr &expr)
		{
			expr_ = &expr;
			return TRAVERSE_STOP;
		}

		FindExpression()
		: expr_(nullptr)
		{}

		inline bool found() const
		{
			return expr_ != nullptr;
		}

		inline const Expr &get() const
		{
			return *expr_;
		}

		template<class ExprTree>
		inline bool apply(const ExprTree &expr)
		{
			traverse(expr, *this);
			return found();
		}

		const Expr * expr_;
	};


	template<template<class...> class Expr>
	class TPLExpressionExists {
	public:

		template<class Any>
		inline constexpr static int visit(const Any &) { return TRAVERSE_CONTINUE; }
		
		template<class... Inner>
		int visit(const Expr<Inner...> &expr)
		{
			found_ = true;
			return TRAVERSE_STOP;
		}

		TPLExpressionExists()
		: found_(false)
		{}

		inline bool found() const
		{
			return found_;
		}


		template<class ExprTree>
		inline bool apply(const ExprTree &expr)
		{
			traverse(expr, *this);
			return found();
		}

		bool found_;
	};

	template<class FunctionSpaceT, class ExprTree>
	inline std::shared_ptr<FunctionSpaceT> trial_space(const ExprTree &tree)
	{
		std::shared_ptr<FunctionSpaceT> ret = nullptr;

		FindExpression< TrialFunction<FunctionSpaceT> > f;
		if(f.apply(tree)) {
			ret = f.get().space_ptr();
		} else {
			std::cerr << "[Warning] unable to find TrialFunction in " << (tree.getClass()) << std::endl; 
		}


		return ret;
	}

	template<class FunctionSpaceT, class ExprTree>
	inline std::shared_ptr<FunctionSpaceT> test_space(const ExprTree &tree)
	{
		std::shared_ptr<FunctionSpaceT> ret = nullptr;

		FindExpression< TestFunction<FunctionSpaceT> > f;
		if(f.apply(tree)) {
			ret = f.get().space_ptr();
		} else {
			std::cerr << "[Warning] unable to find TestFunction in " << (tree.getClass()) << std::endl; 
		}

		return ret;
	}

	template<class ExprTree>
	inline bool is_trial(const ExprTree &tree)
	{
		TPLExpressionExists<utopia::TrialFunction> f;
		return f.apply(tree);
	}

	template<class ExprTree>
	inline bool is_test(const ExprTree &tree)
	{
		TPLExpressionExists<utopia::TestFunction> f;
		return f.apply(tree);
	}
}

#endif //UTOPIA_TRAVERSE_HPP

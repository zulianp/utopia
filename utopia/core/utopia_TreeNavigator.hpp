#ifndef UTOPIA_TREE_NAVIGATOR_HPP
#define UTOPIA_TREE_NAVIGATOR_HPP 

#include "utopia_Expressions.hpp"

// #ifdef WITH_OPENCL
// #include "utopia_Evaluate.hpp"
// #endif //WITH_OPENCL


namespace utopia {
	template<class Action>
	class TreeNavigator {
	public:
		template<class Derived>
		void visit(const Expression<Derived> &expr)
		{
		#ifndef NDEBUG				
			ScopedRecursionCounter src(n_recursions_);
			assert(src.good() && "TreeNavigator: accept method must be implented in Derived class");
		#endif //NDEBUG
			visit(expr.derived());
		}

		template<class InnerExpr, class Operation>
		void visit(const Unary<InnerExpr, Operation> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);
			
			visit(expr.expr());
			
			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}

		template<class InnerExpr>
		void visit(const Transposed<InnerExpr> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);
			
			visit(expr.expr());
			
			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}

		template<typename T>
		void visit(const Number<T> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);
			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}

		template<class Left, class Right, class Operation>
		void visit(const Binary<Left, Right, Operation> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);
			
			visit(expr.left());

			action_.in_order_visit(expr);

			visit(expr.right());
			
			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}

		template<class Left, class Right>
		void visit(const Multiply<Left, Right> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);
			
			visit(expr.left());

			action_.in_order_visit(expr);

			visit(expr.right());
			
			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}

		template<class Tensor, int Order>
		void visit(const Wrapper<Tensor, Order> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);
			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}

		template<class Expr, int Number>
		void visit(const Variable<Expr, Number> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);
			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}

// #ifdef WITH_OPENCL
		template<class Expr, int Order>
		void visit(const Evaluate<Expr, Order> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);
			
			//evaluate is a non terminal symbol
			if(!prune_evaluate_branch_) visit(expr.expr());

			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}
// #endif //WITH_OPENCL		

		template<class InnerExpr, class Operation>
		void visit(const Reduce<InnerExpr, Operation> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);

			visit(expr.expr());
			
			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}

		template<class Left, class Right>
		void visit(const Construct<Left, Right> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);
			
			visit(expr.left());
			
			action_.in_order_visit(expr);

			visit(expr.right());
			
			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}

		template<class Type, int Order>
		void visit(const Factory<Type, Order> &expr)
		{
			pre_intercept(expr);

			action_.pre_order_visit(expr);
			action_.post_order_visit(expr);
			
			post_intercept(expr);
		}

		class ScopedRecursionCounter {
		public:
			inline ScopedRecursionCounter(int &n_recursions)
			: n_recursions_(n_recursions)
			{
				n_recursions_;
			}

			inline ~ScopedRecursionCounter()
			{
				--n_recursions_;
			}

			inline bool good() const
			{
				return n_recursions_ == 1;
			}

		private:
			int &n_recursions_;
		};

		TreeNavigator(Action action)
		: action_(action), n_recursions_(0), verbose_(false), prune_evaluate_branch_(false)
		{}

		void setVerbose(const bool verbose)
		{
			verbose_ = verbose;
		}

		void set_prune_evaluate_branch(const bool value)
		{
			prune_evaluate_branch_ = value;
		}

		template<class Expr>
		void pre_intercept(const Expr &expr)
		{

		}

		template<class Expr>
		void post_intercept(const Expr &expr)
		{
			if(verbose_) std::cout << "visited " << expr.getClass() << std::endl;
		}

	private:
		Action action_;
		int n_recursions_;
		bool verbose_;
		bool prune_evaluate_branch_;
	};

	template<class Action>
	TreeNavigator<const Action &> make_nav(const Action &action)
	{
		return TreeNavigator<const Action &>(action);
	}

	template<class Action>
	TreeNavigator<Action &> make_nav(Action &action)
	{
		return TreeNavigator<Action &>(action);
	}
}

#endif //UTOPIA_TREE_NAVIGATOR_HPP

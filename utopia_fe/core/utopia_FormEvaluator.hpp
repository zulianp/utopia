#ifndef UTOPIA_FE_FORM_EVALUATOR_HPP
#define UTOPIA_FE_FORM_EVALUATOR_HPP

#include "utopia_FormExpressions.hpp"
#include "utopia_Traits.hpp"
#include "utopia_fe_lang.hpp"

namespace utopia {

	template<int BAKEND_FLAG>
	class FormContext {
	public:
		template<class Expr>
		void init(const Expr &) {}

		template<class Tensor>
		void init_tensor(Tensor &, const bool) {}
	};

	template<class Form, int BAKEND_FLAG>
	class FormEval {
	public:
		FormEval()
		{
			static_assert(BAKEND_FLAG < utopia::HOMEMADE, "FormEval: unimplemented evaluator for backend");
		}

		template<class Expr, class Tensor, class Context>
		static void apply(const Expr &expr, Tensor &, const Context &) {
			std::cerr << "[Error] not implemented" << std::endl;
			std::cerr << tree_format(expr.getClass()) << std::endl;
		}
	};

	template<class Form, int BAKEND_FLAG = Traits<Form>::Backend>
	class FormEvaluator {
	public:
		template<class Tensor, int Order>
		static void eval(const Integral<Form> &expr, Wrapper<Tensor, Order> &tensor, const bool reset_tensor)
		{	
			FormContext<BAKEND_FLAG> ctx;
			ctx.init(expr);
			ctx.init_tensor(tensor, reset_tensor);
			FormEval<Form, BAKEND_FLAG>::apply(expr, tensor, ctx);
		}

		template<class Tensor, int Order>
		static void eval(const Integral<Form> &expr, Wrapper<Tensor, Order> &tensor, FormContext<BAKEND_FLAG> &ctx, const bool reset_tensor)
		{	
			ctx.init_tensor(tensor, reset_tensor);
			FormEval<Form, BAKEND_FLAG>::apply(expr, tensor, ctx);
		}
	};
}

#endif //UTOPIA_FE_FORM_EVALUATOR_HPP

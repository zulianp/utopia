#ifndef UTOPIA_FE_FUNCTION_HPP
#define UTOPIA_FE_FUNCTION_HPP 

#include <functional>

namespace utopia {

	template<class Fun>
	class ContextFunction : public Expression< ContextFunction<Fun> > {
	public:

		typedef utopia::Traits<typename Fun::result_type> Traits;
		typedef typename Traits::Scalar Scalar;
		static const int Order = Traits::Order;

		template<int Backend>
		auto eval(const AssemblyContext<Backend> &ctx) -> decltype(fun_(ctx))
		{
			return fun_(ctx);
		}

		ContextFunction(Fun fun)
		: fun_(fun)
		{}

		std::string getClass() const { return "ContextFunction"; }

	private:
		Fun fun_;
	};


	template<class Fun>
	class Traits< ContextFunction<Fun> > : public Traits<typename Fun::result_type> {};

	template<class Fun>
	inline ContextFunction<Fun> ctx_fun(Fun f)
	{
		return ContextFunction<Fun>(f);
	}
}

#endif //UTOPIA_FE_FUNCTION_HPP

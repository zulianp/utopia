#ifndef UTOPIA_FE_EXPRESSION_HPP
#define UTOPIA_FE_EXPRESSION_HPP 

namespace utopia {
	class FEExpression {
	public:
		virtual ~FEExpression() {}
		inline static bool is_fe() { return true; }
	};
}

#endif //UTOPIA_FE_EXPRESSION_HPP

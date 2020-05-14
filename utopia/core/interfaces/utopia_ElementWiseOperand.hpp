#ifndef UTOPIA_ELEMENT_WISE_OPERAND_HPP
#define UTOPIA_ELEMENT_WISE_OPERAND_HPP

namespace utopia {
	template<class T>
	class ElementWiseOperand {
	public:
            virtual ~ElementWiseOperand() = default;
            virtual void e_mul(const T &other) = 0;
            virtual void e_div(const T &other) = 0;
            virtual void e_min(const T &other) = 0;
            virtual void e_max(const T &other) = 0;
	};
}

#endif //UTOPIA_ELEMENT_WISE_OPERAND_HPP
#ifndef UTOPIA_CONDITIONAL_TYPE_HPP
#define UTOPIA_CONDITIONAL_TYPE_HPP

namespace utopia {
	class NullType {};

	template<class Type, int True>
	class ConditionalType {};


	template<class Type_>
	class ConditionalType<Type_, 1> {
	public:
		using Type = Type_;
	};

	template<class Type_>
	class ConditionalType<Type_, 0> {
	public:
		using Type = NullType;
	};

}

#endif //UTOPIA_CONDITIONAL_TYPE_HPP

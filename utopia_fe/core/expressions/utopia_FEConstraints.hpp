#ifndef UTOPIA_FE_CONSTARAINTS_HPP
#define UTOPIA_FE_CONSTARAINTS_HPP 

#include "utopia_Equations.hpp"

namespace utopia {
	template<class... Constraint>
	class FEConstraints {
	public:

		static const int n_constraints = std::tuple_size< std::tuple<Constraint...> >::value;

		FEConstraints(const Constraint &...eqs)
		: eqs_(eqs...)
		{ }

		FEConstraints(const std::tuple<Constraint...> &eqs)
		: eqs_(eqs)
		{}

		template<int Index>
		inline auto get() const -> const typename std::tuple_element<Index, std::tuple<Constraint...>>::type
		{
			return std::get<Index>(eqs_);
		}

		template<int Index>
		inline auto get() -> typename std::tuple_element<Index, std::tuple<Constraint...>>::type
		{
			return std::get<Index>(eqs_);
		}

		template<class Fun>
		void each(Fun fun) 
		{
			EquationIterator<FEConstraints, 0, n_constraints> iter(*this);
			iter.visit(fun);
		}

		template<class Fun>
		void each(Fun fun) const
		{
			EquationIterator<const FEConstraints, 0, n_constraints> iter(*this);
			iter.visit(fun);
		}

		const std::tuple<Constraint...> &equations() const
		{
			return eqs_;
		}

	private:
		std::tuple<Constraint...> eqs_;
	};

	template<>
	class FEConstraints<> {
	public:
		template<class Fun>
		void each(Fun fun) const
		{

		}
	};

	template<class... Constraint>
	inline FEConstraints<Constraint...> constraints(const Constraint &... eqs)
	{
		return FEConstraints<Constraint...>(eqs...);
	}

	template<class... Constraint, class Appended>
	FEConstraints<Constraint..., Appended> operator+(const FEConstraints<Constraint...> &constr, const Appended &eq)
	{
		return FEConstraints<Constraint..., Appended>(std::tuple_cat(constr.equations(), std::make_tuple(eq)));
	}
}

#endif //UTOPIA_FE_CONSTARAINTS_HPP

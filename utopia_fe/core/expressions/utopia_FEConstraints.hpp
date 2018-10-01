#ifndef UTOPIA_FE_CONSTARAINTS_HPP
#define UTOPIA_FE_CONSTARAINTS_HPP

#include "utopia_Equations.hpp"

#include "tao/tuple/tuple.hpp"

namespace utopia {
	template<class... Constraint>
	class FEConstraints {
	public:

		static const int n_constraints = sizeof...(Constraint);

		FEConstraints(const Constraint &...eqs)
		: eqs_(eqs...)
		{ }

		FEConstraints(const tao::tuple<Constraint...> &eqs)
		: eqs_(eqs)
		{}

		template<int Index>
		inline auto get() const -> const typename tao::tuple_element<Index, tao::tuple<Constraint...>>::type
		{
			return tao::get<Index>(eqs_);
		}

		template<int Index>
		inline auto get() -> typename tao::tuple_element<Index, tao::tuple<Constraint...>>::type
		{
			return tao::get<Index>(eqs_);
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

		const tao::tuple<Constraint...> &equations() const
		{
			return eqs_;
		}

	private:
		tao::tuple<Constraint...> eqs_;
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
		return FEConstraints<Constraint..., Appended>(tao::tuple_cat(constr.equations(), tao::make_tuple(eq)));
	}
}

#endif //UTOPIA_FE_CONSTARAINTS_HPP

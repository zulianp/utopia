#ifndef UTOPIA_EQUATIONS_HPP
#define UTOPIA_EQUATIONS_HPP 

#include <utility>

namespace utopia {

	template<class Eqs, int Begin, int End>
	class EquationIterator {
	public:
		template<class Fun>
		void visit(Fun fun)
		{
			fun(Begin, eqs.template get<Begin>());
			EquationIterator<Eqs, Begin + 1, End> next(eqs);
			next.visit(fun);
		}

		EquationIterator(Eqs &eqs) : eqs(eqs) {}
		Eqs &eqs;
	};

	template<class Eqs, int Begin>
	class EquationIterator<Eqs, Begin, Begin + 1> {
	public:
		template<class Fun>
		void visit(Fun fun)
		{
			fun(Begin, eqs.template get<Begin>());
		}

		EquationIterator(Eqs &eqs) : eqs(eqs) {}
		Eqs &eqs;
	};

	template<class... Equation>
	class Equations {
	public:

		static const int n_equations = std::tuple_size< std::tuple<Equation...> >::value;

		Equations(const Equation &...eqs)
		: eqs_(eqs...)
		{ }

		template<int Index>
		inline auto get() const -> const typename std::tuple_element<Index, std::tuple<Equation...>>::type
		{
			return std::get<Index>(eqs_);
		}

		template<int Index>
		inline auto get() -> typename std::tuple_element<Index, std::tuple<Equation...>>::type
		{
			return std::get<Index>(eqs_);
		}

		template<class Fun>
		void each(Fun fun)
		{
			EquationIterator<Equations, 0, n_equations> iter(*this);
			iter.visit(fun);
		}

		template<class Fun>
		void each(Fun fun) const
		{
			EquationIterator<const Equations, 0, n_equations> iter(*this);
			iter.visit(fun);
		}

	private:
		std::tuple<Equation...> eqs_;
	};


	template<class... Equation>
	inline Equations<Equation...> equations(const Equation &... eqs)
	{
		return Equations<Equation...>(eqs...);
	}

	// //base case
	// template<class Left, class Right>
	// class Traits<Equations<Left, Right>> {
	// public:
	// 	typedef typename utopia::FormTraits<Left>::Implementation Implementation;
	// 	typedef typename utopia::FormTraits<Left>::Scalar Scalar;

	// 	static const int Order   = utopia::FormTraits<Left>::Order + FormTraits<Right>::Order;
	// 	static const int Backend = utopia::FormTraits<Left>::Backend;
	// };

	// //variadic case
	// template<class... Equation>
	// class Traits< Equations<Equation...> > {
	// public:
	// 	template<class Eqs>
	// 	inline static constexpr int order()
	// 	{
	// 		return Traits<Eqs>::Order;
	// 	}

	// 	template<class First, class...Rest>
	// 	inline static constexpr int order()
	// 	{
	// 		return Traits<First>::Order + order<Rest...>();
	// 	}

		template<class First, class...Rest>
		class GetFirst {
		public:
			typedef First Type;
		};

	// 	typedef typename GetFirst<Equation...>::Type First;
	// 	typedef typename utopia::FormTraits<First>::Implementation Implementation;
	// 	typedef typename utopia::FormTraits<First>::Scalar Scalar;

	// 	const static int Order = order<Equation...>();
	// 	static const int Backend = utopia::FormTraits<First>::Backend;
	// };
}

#endif //UTOPIA_EQUATIONS_HPP

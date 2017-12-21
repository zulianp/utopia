#ifndef UTOPIA_OPTIONAL_HPP
#define UTOPIA_OPTIONAL_HPP 

#include <utility>

namespace utopia {

	template<class Args, int Begin, int End>
	class OptionalIterator {
	public:
		template<class Fun>
		void visit(Fun &fun)
		{
			fun.parse_arg(args.template get<Begin>());
			OptionalIterator<Args, Begin + 1, End> next(args);
			next.visit(fun);
		}

		OptionalIterator(Args &args) : args(args) {}
		Args &args;
	};

	template<class Args, int Begin>
	class OptionalIterator<Args, Begin, Begin> {
	public:
		template<class Fun>
		void visit(Fun &) { }
		OptionalIterator(const Args &){}
	};

	//concept
	template<class... Args>
	class Optional {
	public:
		static const int n_args = std::tuple_size< std::tuple<Args...> >::value;

		Optional(const Args &...args)
		: opts(std::make_tuple(args...))
		{}

		template<int Index>
		inline auto get() -> typename std::tuple_element<Index, std::tuple<Args...>>::type
		{
			return std::get<Index>(opts);
		}

		template<int Index>
		inline auto get() const -> const typename std::tuple_element<Index, std::tuple<Args...>>::type
		{
			return std::get<Index>(opts);
		}


		template<class Fun>
		void each(Fun &fun)
		{
			OptionalIterator<Optional, 0, n_args> iter(*this);
			iter.visit(fun);
		}

		template<class Fun>
		void each(Fun &fun) const
		{
			OptionalIterator<const Optional, 0, n_args> iter(*this);
			iter.visit(fun);
		}

		std::tuple<Args...> opts;
	};

	template<class... Args>
	Optional<typename std::remove_reference<Args>::type...> options(Args &&...args)
	{
		return Optional<typename std::remove_reference<Args>::type...>(std::forward<Args...>(args...));
	} 
}

#endif //UTOPIA_OPTIONAL_HPP

// #ifndef UTOPIA_MIXED_FE_SPACE_HPP
// #define UTOPIA_MIXED_FE_SPACE_HPP

// #include <iostream>
// #include <tuple>

// namespace utopia {

// 	template<class Tuple, int Index, int N>
// 	struct TupleEach
// 	{
// 		template<class Functor>
// 		static void apply(Tuple &tuple, Functor &&fun) {

// 			fun. template call<Index>(std::get<Index>(tuple));

// 			TupleEach<Tuple, Index+1, N>::apply(tuple, std::forward<Functor>(fun));
// 		}
// 	};

// 	template<class Tuple, int N>
// 	struct TupleEach<Tuple, N, N>
// 	{
// 		template<class Functor>
// 		static void apply(Tuple &, Functor &&) {}
// 	};

// 	//example functor
// 	struct PrintIndex {
// 		template<int Index, class Arg>
// 		void call(Arg &arg) const
// 		{
// 			std::cout << "Index:" << Index << std::endl;
// 		}
// 	};

// 	template<class ...Args>
// 	class MixedFESpace {
// 	public:
// 		typedef std::tuple<Args...> Spaces;
// 		typedef std::tuple<typename Args::Function...> Functions;

// 		MixedFESpace(Args &...args)
// 		: spaces_(args...)
// 		{}

// 		MixedFESpace(Spaces &&spaces)
// 		: spaces_(spaces)
// 		{}

// 		template<class Functor>
// 		void each(Functor &&fun)
// 		{
// 			TupleEach<Spaces, 0, std::tuple_size<Spaces>::value>::apply(spaces_, fun);
// 		}

// 		inline Spaces &subspaces()
// 		{
// 			return spaces_;
// 		}
// 		inline static constexpr std::size_t n_subspaces()
// 		{
// 			return std::tuple_size<Spaces>::value;
// 		}

// 		template<std::size_t Index>
// 		typename std::tuple_element<Index, Spaces>::type &sub()
// 		{
// 			return std::get<Index>(spaces_);
// 		}

// 	private:
// 		Spaces spaces_;
// 	};

// 	template<class First, class... Rest>
// 	class MixedFEFunction : public FEFunction<typename First::Scalar> {
// 	public:
// 		MixedFEFunction(MixedFESpace<First, Rest...> &space)
// 		{}
// 	};

// 	template<class... Args>
// 	MixedFEFunction<Args...> fe_function(MixedFESpace<Args...> &space)
// 	{
// 		return MixedFEFunction<Args...>(space);
// 	}

// 	template<class Traits1, class Traits2, int Backend>
// 	MixedFESpace<FESpace<Traits1, Backend>, FESpace<Traits2, Backend> > operator*(FESpace<Traits1, Backend> &fs1,
// 																				  FESpace<Traits2, Backend> &fs2)

// 	{
// 		return MixedFESpace<FESpace<Traits1, Backend>, FESpace<Traits2, Backend> >(fs1, fs2);
// 	}

// 	template<class... Spaces, class Traits, int Backend>
// 	MixedFESpace<Spaces..., FESpace<Traits, Backend> > operator*(MixedFESpace<Spaces...>  &&fs1,
// 																 FESpace<Traits, Backend> &fs2)

// 	{
// 		return MixedFESpace<Spaces..., FESpace<Traits, Backend> >(std::tuple_cat(fs1.subspaces(), std::tie(fs2)));
// 	}

// 	template<class... Spaces, class Traits, int Backend>
// 	MixedFESpace<FESpace<Traits, Backend>, Spaces...> operator*(FESpace<Traits, Backend> &fs1,
// 																MixedFESpace<Spaces...>  &&fs2)

// 	{
// 		return MixedFESpace<FESpace<Traits, Backend>, Spaces...>(std::tuple_cat(std::tie(fs1), fs2.subspaces()));
// 	}

// 	template<class... Spaces1, class... Spaces2>
// 	MixedFESpace<Spaces1..., Spaces2...> operator*(MixedFESpace<Spaces1...>  &&fs1,
// 												   MixedFESpace<Spaces2...>  &&fs2)

// 	{
// 		return MixedFESpace<Spaces1..., Spaces2...>( std::tuple_cat(fs1.subspaces(), fs2.subspaces()) );
// 	}

// }

// #endif //UTOPIA_MIXED_FE_SPACE_HPP


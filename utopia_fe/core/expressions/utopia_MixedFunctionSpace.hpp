#ifndef UTOPIA_MIXED_FUNCTION_SPACE_HPP
#define UTOPIA_MIXED_FUNCTION_SPACE_HPP

#include <utility>
#include "tao/tuple/tuple.hpp"

namespace utopia {

	template<class InputTuple,
			 class OutputTuple,
			 std::size_t N>
	struct TupleFilter {

		template<class Function>
	    static void apply(const InputTuple& in, OutputTuple& out, Function fun)
	    {
	        TupleFilter<InputTuple, OutputTuple, N-1>::apply(in, out, fun);
	        fun(tao::get<N-1>(in), tao::get<N-1>(out));
	    }
	};

	template<class InputTuple,
			 class OutputTuple>
	struct TupleFilter<InputTuple, OutputTuple, 1> {

		template<class Function>
	    static void apply(const InputTuple& in, OutputTuple& out, Function fun)
	    {
	        fun(tao::get<0>(in), tao::get<0>(out));
	    }
	};

	template<class...T>
	using MixedFunctionSpace = tao::tuple<T...>;


	class MakeTrial {
	public:
		template<class Space, class Trial>
		void operator()(const Space &s, Trial &t)
		{
			s = trial(t);
		}
	};

	class MakeTest {
	public:
		template<class Space, class Test>
		void operator()(const Space &s, Test &t)
		{
			s = test(t);
		}
	};

	template<class... Args>
	MixedFunctionSpace<FunctionSpace<Args>...> mixed(FunctionSpace<Args> &&...args)
	{
		return tao::make_tuple(args...);
	}

	template<class... Args>
	MixedFunctionSpace<std::shared_ptr<FunctionSpace<Args>>...> mixed(std::shared_ptr<FunctionSpace<Args> > &&... args)
	{
		return tao::make_tuple(args...);
	}

	template<class... Args>
	MixedFunctionSpace<TrialFunction<FunctionSpace<Args>>...> trials(MixedFunctionSpace<FunctionSpace<Args>...> &&mixed_space)
	{
		MixedFunctionSpace<TrialFunction<FunctionSpace<Args>>...> ret;
		TupleFilter<decltype(mixed_space), decltype(ret), decltype(mixed_space)::tuple_size>::apply(
			mixed_space,
			ret,
			MakeTrial()
		);

		return ret;
	}

	template<class... Args>
	MixedFunctionSpace<TestFunction<FunctionSpace<Args>>...> tests(MixedFunctionSpace<FunctionSpace<Args>...> &&mixed_space)
	{
		MixedFunctionSpace<TestFunction<FunctionSpace<Args>>...> ret;
		TupleFilter<decltype(mixed_space), decltype(ret), decltype(mixed_space)::tuple_size>::apply(
			mixed_space,
			ret,
			MakeTest()
		);

		return ret;
	}

}

#endif //UTOPIA_MIXED_FUNCTION_SPACE_HPP

// #ifndef UTOPIA_TRILINOS_EVAL_DISTANCE_HPP
// #define UTOPIA_TRILINOS_EVAL_DISTANCE_HPP

// #include "utopia_Eval_Empty.hpp"
// #include "utopia_kokkos_Eval_Distance.hpp"

// namespace utopia {

// 	namespace trilinos_ {
// 		template<class Left, class Right, int NormType, class Traits>
// 		class EvalDistance {
// 		public:
// 		    typedef typename Traits::Scalar Scalar;
// 		    typedef typename Traits::Vector Vector;

// 		    inline static Scalar apply(const Left &left, const Right &right)
// 		    {
// 		        Scalar result;

// 		        auto &&l = Eval<Left, Traits>::apply(left);
// 		        auto &&r = Eval<Right, Traits>::apply(right);

// 		        result = KokkosEvalDistance<Vector, NormType>::apply(l, r, false);

// 		        auto &comm = *l.communicator();
// 		        Scalar result_global = 0.;
// 		        Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &result, &result_global);

// 		        return KokkosEvalDistance<Vector, NormType>::finalize(result_global);
// 		    }
// 		};
// 	}

// 	template<class Tensor_, class Traits>
// 	class Eval<
// 		Distance<
// 			Wrapper<Tensor_, 1>,
// 			Wrapper<Tensor_, 1>,
// 			2
// 			>,
// 			Traits, TRILINOS> {
// 	public:
// 		using Tensor = utopia::Wrapper<Tensor_, 1>;
// 		using Expr = utopia::Distance<Tensor, Tensor, 2>;

// 	    typedef typename Traits::Scalar Scalar;

// 	    inline static Scalar apply(const Expr &expr)
// 	    {
// 	        UTOPIA_TRACE_BEGIN(expr);

// 	     	Scalar result = trilinos_::EvalDistance<Tensor, Tensor, 2, Traits>::apply(
// 	     		expr.expr().left(),
// 	     		expr.expr().right()
// 	     		);

// 	        UTOPIA_TRACE_END(expr);
// 	        return result;
// 	    }
// 	};

// 	template<class Tensor_, class Traits>
// 	class Eval<
// 		Distance<
// 			Wrapper<Tensor_, 1>,
// 			Wrapper<Tensor_, 1>,
// 			1
// 			>,
// 			Traits, TRILINOS> {
// 	public:
// 		using Tensor = utopia::Wrapper<Tensor_, 1>;
// 		using Expr = utopia::Distance<Tensor, Tensor, 1>;

// 	    typedef typename Traits::Scalar Scalar;

// 	    inline static Scalar apply(const Expr &expr)
// 	    {
// 	        Scalar result;
// 	        UTOPIA_TRACE_BEGIN(expr);

// 	     	result = trilinos_::EvalDistance<Tensor, Tensor, 1, Traits>::apply(
// 	     		expr.expr().left(),
// 	     		expr.expr().right()
// 	     	);

// 	        UTOPIA_TRACE_END(expr);
// 	        return result;
// 	    }
// 	};

// 	template<class Tensor_, class Traits>
// 	class Eval<
// 		Distance<
// 			Wrapper<Tensor_, 1>,
// 			Wrapper<Tensor_, 1>,
// 			INFINITY_NORM_TAG
// 			>,
// 			Traits, TRILINOS> {
// 	public:
// 		using Tensor = utopia::Wrapper<Tensor_, 1>;
// 		using Expr = utopia::Distance<Tensor, Tensor, 1>;

// 	    typedef typename Traits::Scalar Scalar;

// 	    inline static Scalar apply(const Expr &expr)
// 	    {
// 	        Scalar result;
// 	        UTOPIA_TRACE_BEGIN(expr);

// 	     	result = trilinos_::EvalDistance<Tensor, Tensor, INFINITY_NORM_TAG, Traits>::apply(
// 	     		expr.expr().left(),
// 	     		expr.expr().right()
// 	     	);

// 	        UTOPIA_TRACE_END(expr);
// 	        return result;
// 	    }
// 	};

// }

// #endif //UTOPIA_TRILINOS_EVAL_DISTANCE_HPP

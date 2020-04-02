#ifndef UTOPIA_MAT_CHOP_HPP
#define UTOPIA_MAT_CHOP_HPP

#include "utopia_Traits.hpp"
#include "utopia_Algorithms.hpp"

namespace utopia {

	template<class Matrix, int Backend = Traits<Matrix>::Backend>
	class Chop {
	public:
	    using Scalar   = typename utopia::Traits<Matrix>::Scalar;
	    using SizeType = typename utopia::Traits<Matrix>::SizeType;

	    static void apply(Matrix &mat, const Scalar &eps)
	    {
	        each_transform(mat, [eps](const SizeType &, const SizeType &, const Scalar &v) -> Scalar {
	            return device::abs(v) < eps ? 0.0 : v;
	        });
	    }
	};

	template<class Matrix>
	void chop(Tensor<Matrix, 2> &A, const typename Traits<Matrix>::Scalar &eps)
	{
	    Chop<Matrix>::apply(A.derived(), eps);
	}

	template<class Matrix, int Backend = Traits<Matrix>::Backend>
	class ChopSmallerThan {
	public:
	    using Scalar   = typename utopia::Traits<Matrix>::Scalar;
	    using SizeType = typename utopia::Traits<Matrix>::SizeType;

	    static void apply(Matrix &mat, const Scalar &eps)
	    {
	        each_transform(mat, [eps](const SizeType &, const SizeType &, const Scalar &v) -> Scalar {
	            return v < eps ? 0.0 : v;
	        });
	    }
	};

	template<class Matrix, int Backend = Traits<Matrix>::Backend>
	class ChopGreaterThan {
	public:
	    using Scalar   = typename utopia::Traits<Matrix>::Scalar;
	    using SizeType = typename utopia::Traits<Matrix>::SizeType;

	    static void apply(Matrix &mat, const Scalar &eps)
	    {
	        each_transform(mat, [eps](const SizeType &, const SizeType &, const Scalar &v) -> Scalar {
	            return v > eps ? 0.0 : v;
	        });
	    }
	};

	template<class Matrix>
	void chop_smaller_than(Tensor<Matrix, 2> &A, const double eps)
	{
	    ChopSmallerThan<Matrix>::apply(A.derived(), eps);
	}

	template<class Matrix>
	void chop_greater_than(Tensor<Matrix, 2> &A, const double eps)
	{
	    ChopGreaterThan<Matrix>::apply(A.derived(), eps);
	}

}


#endif //UTOPIA_MAT_CHOP_HPP

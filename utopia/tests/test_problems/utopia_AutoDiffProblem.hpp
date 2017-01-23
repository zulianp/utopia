#ifndef UTOPIA_AUTO_DIFF_PROBLEM_HPP
#define UTOPIA_AUTO_DIFF_PROBLEM_HPP 


#include "utopia_Function.hpp"
#include "utopia_AutoDiff.hpp"

namespace utopia {

	template<class Matrix, class Vector, class Fun>
	class AutoDiffProblem : public Function<Matrix, Vector> {
	public:
	    DEF_UTOPIA_SCALAR(Matrix);


	    inline bool value(const Vector &point, Scalar &result) const override {
	        result = scalar_cast<Scalar>(fun_(point));
	        return true;
	    }

	    inline bool gradient(const Vector &point, Vector &result) const override {
	        result = derivative(fun_(independent_variable(point)));
	        return true;
	    }

	    inline bool hessian(const Vector &point, Matrix &result) const override {
	        result = derivative(derivative(fun_(independent_variable(point))));
	        return true;
	    }

	    AutoDiffProblem(Fun fun) : fun_(fun) { }

	private:
		Fun fun_;
	};

	template<class Matrix, class Vector, class Fun>
	AutoDiffProblem<Matrix, Vector, Fun> auto_diff_fun(Fun fun)
	{
		return AutoDiffProblem<Matrix, Vector, Fun>(fun);
	}
}

#endif //UTOPIA_AUTO_DIFF_PROBLEM_HPP

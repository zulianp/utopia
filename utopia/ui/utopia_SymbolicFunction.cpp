#include "utopia_SymbolicFunction.hpp"

#ifdef WITH_TINY_EXPR

#include "tinyexpr.h"
#include "utopia_make_unique.hpp"

#include <vector>

namespace utopia {
	class SymbolicFunction::Impl {
	public:
		Impl(const std::string &expr)
		:
		x_(0.),
		y_(0.),
		z_(0.),
		expr_(expr),
		vars_({
			{"x", &x_},
			{"y", &y_},
			{"z", &z_},
		}),
		err_(0)
		{
		 	/* Compile the expression with variables. */
			e_ = te_compile(expr_.c_str(), &vars_[0], vars_.size(), &err_);
		}

		inline bool valid() const
		{
			return err_ == 0;
		}

		~Impl()
		{
			te_free(e_);
		}

		double eval()
		{
			return te_eval(e_);
		}

	public:
		double x_, y_, z_;

	private:
		std::string expr_;

		std::vector<te_variable> vars_;
		int err_;

		te_expr * e_;
	};

	SymbolicFunction::~SymbolicFunction() {}

	SymbolicFunction::SymbolicFunction(const std::string &expr)
	{
		impl_ = make_unique<Impl>(expr);
	}

	bool SymbolicFunction::valid() const
	{
		return impl_->valid();
	}

	double SymbolicFunction::eval(const double x)
	{
		impl_->x_ = x;
		impl_->y_ = 0.;
		impl_->z_ = 0.;
		return  impl_->eval();
	}

	double SymbolicFunction::eval(const double x, const double y)
	{
		impl_->x_ = x;
		impl_->y_ = y;
		impl_->z_ = 0.;
		return  impl_->eval();
	}

	double SymbolicFunction::eval(const double x, const double y, const double z)
	{
		impl_->x_ = x;
		impl_->y_ = y;
		impl_->z_ = z;
		return  impl_->eval();
	}

	double SymbolicFunction::eval(const std::vector<double> &x)
	{
		//FIXME
		assert(x.size() <= 3);

		switch(x.size()) {
			case 1: {
				return eval(x[0]);
			}

			case 2: {
				return eval(x[0], x[1]);
			}

			case 3: {
				return eval(x[0], x[1], x[2]);
			}

			default:
			{
				std::cerr << "TODO" << std::endl;
				assert(false);
				return 0.;
			}
		}
	}
}

#endif //WITH_TINY_EXPR

#ifndef UTOPIA_SYMBOLIC_FUNCTION_HPP
#define UTOPIA_SYMBOLIC_FUNCTION_HPP

#include "utopia_Base.hpp"

#ifdef WITH_TINY_EXPR

#include <memory>

namespace utopia {

	class SymbolicFunction {
	public:
		class Impl;

		double eval(const double x);
		double eval(const double x, const double y);
		double eval(const double x, const double y, const double z);
		double eval(const double x, const double y, const double z, const double t);
		double eval(const std::vector<double> &x);

		bool valid() const;

		SymbolicFunction(const std::string &expr);
		~SymbolicFunction();

	private:
		std::unique_ptr<Impl> impl_;
	};

}

#endif //WITH_TINY_EXPR

#endif //UTOPIA_SYMBOLIC_FUNCTION_HPP

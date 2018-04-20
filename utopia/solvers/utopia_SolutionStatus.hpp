#ifndef UTOPIA_SOLUTION_STATUS_HPP
#define UTOPIA_SOLUTION_STATUS_HPP 

#include <iostream>

namespace utopia {

	class SolutionStatus {
	public:
		int iterates;
		double function_value;
		double gradient_norm;
		double step_norm;
		int reason;

		SolutionStatus()
		: iterates(0),
		  function_value(-1),
		  gradient_norm(-1),
		  step_norm(-1),
		  reason(-1)
		{}

		void clear()
		{
			iterates = 0;
			function_value = -1;
			gradient_norm = -1;
			step_norm = -1;
			reason = -1;
		}

		inline void describe(std::ostream &os) const
		{
			os << "iterates:       " << iterates << "\n";
			os << "function_value: " << function_value << "\n";
			os << "gradient_norm:  " << gradient_norm << "\n";
			os << "step_norm:      " << step_norm << "\n";
			os << "reason:         " << reason << std::endl;
		}
	};
}

#endif //UTOPIA_SOLUTION_STATUS_HPP

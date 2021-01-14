#ifndef UTOPIA_SOLUTION_STATUS_HPP
#define UTOPIA_SOLUTION_STATUS_HPP

#include <iostream>

namespace utopia {

class SolutionStatus {
 public:
  int iterates{0};
  int num_linear_solves{0};
  double function_value{-1};
  double gradient_norm{-1};
  double relative_gradient_norm{-1};
  double step_norm{-1};
  int reason{-1};
  double execution_time{0};

  SolutionStatus()

      = default;

  void clear() {
    iterates = 0;
    num_linear_solves = 0;
    function_value = -1;
    gradient_norm = -1;
    relative_gradient_norm = -1;
    step_norm = -1;
    reason = -1;
    execution_time = 0;
  }

  inline void describe(std::ostream &os) const {
    os << "iterates:       " << iterates << "\n";
    os << "num_linear_solves:       " << num_linear_solves << "\n";
    os << "function_value: " << function_value << "\n";
    os << "gradient_norm:  " << gradient_norm << "\n";
    os << "relative_gradient_norm: " << relative_gradient_norm << " \n";
    os << "execution_time: " << execution_time << " \n";
    os << "step_norm:      " << step_norm << "\n";
    os << "reason:         " << reason << std::endl;
  }
};
}  // namespace utopia

#endif  // UTOPIA_SOLUTION_STATUS_HPP

#ifndef UTOPIA_REDUNDANT_QP_SOLVER_HPP
#define UTOPIA_REDUNDANT_QP_SOLVER_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class RedundantQPSolver {
    public:
        RedundantQPSolver() {
            static_assert(Backend < HOMEMADE, "RedundantLinearSolver not implemented for this backend");
        }
    };
}  // namespace utopia

#endif  // UTOPIA_REDUNDANT_QP_SOLVER_HPP
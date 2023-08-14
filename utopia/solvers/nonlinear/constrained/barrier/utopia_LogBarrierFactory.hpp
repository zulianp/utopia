#ifndef UTOPIA_LOG_BARRIER_FUNCTION_FACTORY_HPP
#define UTOPIA_LOG_BARRIER_FUNCTION_FACTORY_HPP

#include "utopia_BoundedLogBarrierFunction.hpp"
#include "utopia_LogBarrierFunction.hpp"
#include "utopia_ShiftedPenalty.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class PenaltyFactory {
    public:
        ///////////////////////////////////////////////////////
        using LogBarrierFunction = utopia::LogBarrierFunction<Matrix, Vector>;
        using LogBarrierFunctionPtr = std::unique_ptr<LogBarrierFunction>;

        ///////////////////////////////////////////////////////
        // Log barrier functions
        using Penalty = utopia::Penalty<Matrix, Vector>;
        using PenaltyPtr = std::unique_ptr<Penalty>;
        using LogBarrier = utopia::LogBarrier<Matrix, Vector>;
        using BoundedLogBarrier = utopia::BoundedLogBarrier<Matrix, Vector>;
        using ShiftedPenalty = utopia::ShiftedPenalty<Matrix, Vector>;

        static PenaltyPtr new_penalty(const std::string &function_type) {
            PenaltyPtr function;

            if (function_type == "ShiftedPenalty") {
                function = utopia::make_unique<ShiftedPenalty>();
            } else if (function_type == "BoundedLogBarrier") {
                function = utopia::make_unique<BoundedLogBarrier>();
            } else {
                function = utopia::make_unique<LogBarrier>();
            }

            return function;
        }

        static LogBarrierFunctionPtr new_log_barrier_function(const std::string &function_type) {
            if (function_type == "BoundedLogBarrier") {
                return utopia::make_unique<LogBarrierFunction>(utopia::make_unique<BoundedLogBarrier>());
            } else {
                return utopia::make_unique<LogBarrierFunction>(utopia::make_unique<LogBarrier>());
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_FUNCTION_FACTORY_HPP

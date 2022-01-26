#ifndef UTOPIA_LOG_BARRIER_FUNCTION_FACTORY_HPP
#define UTOPIA_LOG_BARRIER_FUNCTION_FACTORY_HPP

#include "utopia_BoundedLogBarrierFunction.hpp"
#include "utopia_LogBarrierFunction.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class LogBarrierFactory {
    public:
        ///////////////////////////////////////////////////////
        using LogBarrierFunction = utopia::LogBarrierFunction<Matrix, Vector>;
        using LogBarrierFunctionPtr = std::unique_ptr<LogBarrierFunction>;

        ///////////////////////////////////////////////////////
        // Log barrier functions
        using LogBarrierBase = utopia::LogBarrierBase<Matrix, Vector>;
        using LogBarrierBasePtr = std::unique_ptr<LogBarrierBase>;
        using LogBarrier = utopia::LogBarrier<Matrix, Vector>;
        using BoundedLogBarrier = utopia::BoundedLogBarrier<Matrix, Vector>;

        static LogBarrierBasePtr new_log_barrier(const std::string &function_type) {
            LogBarrierBasePtr function;

            if (function_type == "BoundedLogBarrier") {
                function = utopia::make_unique<BoundedLogBarrier>();
            } else {
                function = utopia::make_unique<LogBarrier>();
            }

            return function;
        }

        static LogBarrierFunctionPtr new_log_barrier_function(const std::string &function_type) {
            return utopia::make_unique<LogBarrierFunction>(new_log_barrier(function_type));
        }
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_FUNCTION_FACTORY_HPP

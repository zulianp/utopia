#ifndef UTOPIA_LOG_BARRIER_FUNCTION_FACTORY_HPP
#define UTOPIA_LOG_BARRIER_FUNCTION_FACTORY_HPP

#include "utopia_BoundedLogBarrierFunction.hpp"
#include "utopia_LogBarrierFunction.hpp"
#include "utopia_LogBarrierFunctionWithSelection.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class LogBarrierFunctionFactory {
    public:
        ///////////////////////////////////////////////////////
        // Log Barrier function abstraction
        using LogBarrierFunctionBase = utopia::LogBarrierFunctionBase<Matrix, Vector>;
        using LogBarrierFunctionBasePtr = std::unique_ptr<LogBarrierFunctionBase>;

        using LogBarrierFunction = utopia::LogBarrierFunction<Matrix, Vector>;
        using LogBarrierFunctionWithSelection = utopia::LogBarrierFunctionWithSelection<Matrix, Vector>;
        // using BoundedLogBarrierFunction = utopia::BoundedLogBarrierFunction<Matrix, Vector>;

        ///////////////////////////////////////////////////////
        // Log barrier functions
        using LogBarrierBase = utopia::LogBarrierBase<Matrix, Vector>;
        using LogBarrierBasePtr = std::unique_ptr<LogBarrierBase>;
        using LogBarrier = utopia::LogBarrier<Matrix, Vector>;
        using LogBarrierWithSelection = utopia::LogBarrierWithSelection<Matrix, Vector>;

        static LogBarrierBasePtr new_log_barrier(const std::string &function_type) {
            LogBarrierBasePtr function;

            if (function_type == "LogBarrierWithSelection") {
                function = utopia::make_unique<LogBarrierWithSelection>();
            }
            // else if (function_type == "BoundedLogBarrier") {
            //     function = utopia::make_unique<BoundedLogBarrier>();
            // }
            else {
                function = utopia::make_unique<LogBarrier>();
            }

            return function;
        }

        static LogBarrierFunctionBasePtr new_log_barrier_function(const std::string &function_type) {
            LogBarrierFunctionBasePtr function;

            if (function_type == "LogBarrierFunctionWithSelection") {
                function = utopia::make_unique<LogBarrierFunctionWithSelection>();
            }
            // else if (function_type == "BoundedLogBarrierFunction") {
            //     function = utopia::make_unique<BoundedLogBarrierFunction>();
            // }
            else {
                function = utopia::make_unique<LogBarrierFunction>();
            }

            return function;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_FUNCTION_FACTORY_HPP

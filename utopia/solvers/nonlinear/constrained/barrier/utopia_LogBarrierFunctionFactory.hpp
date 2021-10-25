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
        using LogBarrierFunction = utopia::LogBarrierFunction<Matrix, Vector>;
        // using BoundedLogBarrierFunction = utopia::BoundedLogBarrierFunction<Matrix, Vector>;
        using LogBarrierFunctionWithSelection = utopia::LogBarrierFunctionWithSelection<Matrix, Vector>;
        using LogBarrierFunctionBase = utopia::LogBarrierFunctionBase<Matrix, Vector>;
        using LogBarrierFunctionBasePtr = std::unique_ptr<LogBarrierFunctionBase>;

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

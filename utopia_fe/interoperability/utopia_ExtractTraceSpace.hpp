#ifndef UTOPIA_EXTRACT_TRACE_SPACE_HPP
#define UTOPIA_EXTRACT_TRACE_SPACE_HPP

namespace utopia {
    template <class FunctionSpace, class TraceSpace>
    class ExtractTraceSpace {};

    template <class FunctionSpace, class TraceSpace>
    inline void extract_trace_space(const FunctionSpace &in, TraceSpace &out) {
        ExtractTraceSpace<FunctionSpace, TraceSpace>::apply(in, out);
    }
}  // namespace utopia

#endif  // UTOPIA_EXTRACT_TRACE_SPACE_HPP
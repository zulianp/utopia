#ifndef UTOPIA_CONVERT_FUNCTION_SPACE_HPP
#define UTOPIA_CONVERT_FUNCTION_SPACE_HPP
namespace utopia {
    template <class FunctionSpaceIn, class FunctionSpaceOut>
    class ConvertFunctionSpace {};

    // FIXME use convert once we have static polymorphism up
    template <class FunctionSpaceIn, class FunctionSpaceOut>
    inline void convert_function_space(const FunctionSpaceIn &in, FunctionSpaceOut &out) {
        ConvertFunctionSpace<FunctionSpaceIn, FunctionSpaceOut>::apply(in, out);
    }
}  // namespace utopia

#endif  // UTOPIA_CONVERT_FUNCTION_SPACE_HPP

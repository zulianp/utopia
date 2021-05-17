#ifndef UTOPIA_TENSOR_FUNCTION_SPACE_HPP
#define UTOPIA_TENSOR_FUNCTION_SPACE_HPP

#include "utopia_FEForwardDeclarations.hpp"

namespace utopia {
    template <class Space>
    class TensorFunctionSpace : public FunctionSpace<TensorFunctionSpace<Space> > {
    public:
    };
}  // namespace utopia

#endif  // UTOPIA_TENSOR_FUNCTION_SPACE_HPP

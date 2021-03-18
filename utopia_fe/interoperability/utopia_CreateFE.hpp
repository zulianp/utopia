#ifndef UTOPIA_CREATE_FE_HPP
#define UTOPIA_CREATE_FE_HPP

namespace utopia {
    template <class FunctionSpace, class FE>
    class CreateFE {};

    template <class FunctionSpace, class FE>
    inline void create_fe(const FunctionSpace &space, FE &fe, int order = 0) {
        CreateFE<FunctionSpace, FE>::apply(space, fe, order);
    }
}  // namespace utopia

#endif  // UTOPIA_CREATE_FE_HPP
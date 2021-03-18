#ifndef UTOPIA_CREATE_FE_HPP
#define UTOPIA_CREATE_FE_HPP

namespace utopia {
    template <class FunctionSpace, class FE>
    class CreateFE {
    public:
        // apply
    };

    template <class FunctionSpace, class FE>
    inline void create_fe(const FunctionSpace &space, FE &fe) {
        CreateFE<FunctionSpace, FE>::apply(space, fe);
    }
}  // namespace utopia

#endif  // UTOPIA_CREATE_FE_HPP
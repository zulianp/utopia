#ifndef UTOPIA_MARS_FORWARD_DECLARATIONS_HPP
#define UTOPIA_MARS_FORWARD_DECLARATIONS_HPP

namespace utopia {
    namespace mars {
        class Mesh;
        class FunctionSpace;
        class FEAssembler;
        class OmniAssembler;
        class Factory;
        template <class DMesh, typename...>
        class ConcreteFEAssembler;
        template <class DMesh, int Degree_>
        class FEHandler;
    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_FORWARD_DECLARATIONS_HPP

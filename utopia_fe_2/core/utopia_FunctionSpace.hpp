#ifndef UTOPIA_FUNCTION_SPACE_HPP
#define UTOPIA_FUNCTION_SPACE_HPP

#include "utopia_Mesh.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Mesh>
    class IFunctionSpace : public Configurable, public Describable {
    public:
        using Traits = utopia::Traits<Mesh>;
        using Vector = typename Traits::Vector;
        using Matrix = typename Traits::Matrix;

        ~IFunctionSpace() override = default;
        virtual bool write(const Path &path, const Vector &x) const = 0;
        virtual void create_vector(Vector &x) const = 0;
    };

    template <class... T>
    class FunctionSpace {};
}  // namespace utopia

#endif  // UTOPIA_FUNCTION_SPACE_HPP

#ifndef UTOPIA_MESH_HPP
#define UTOPIA_MESH_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"

namespace utopia {

    class IMesh : public Configurable, public Describable {
    public:
        ~IMesh() override = default;
        virtual bool write(const Path &) const { return false; }
    };

    template <class... T>
    class Mesh {};
}  // namespace utopia

#endif  // UTOPIA_MESH_HPP
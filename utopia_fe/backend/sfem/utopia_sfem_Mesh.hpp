#ifndef UTOPIA_SFEM_MESH_HPP
#define UTOPIA_SFEM_MESH_HPP

#include <memory>

#include "utopia_fe_base.hpp"

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

#include "utopia_sfem_ForwardDeclarations.hpp"

namespace utopia {
    template <>
    class Traits<utopia::sfem::Mesh> : public Traits<UVector> {
        // class Traits<utopia::sfem::Mesh> : public Traits<utopia::TpetraMatrix> {
    public:
        using Super = Traits<UVector>;
        // using SideSet = utopia::sfem::SideSet;
        // using AABB = utopia::AABB<std::vector<Super::Scalar>>;
    };

    namespace sfem {
        class Mesh final : public Configurable, public Describable {
        public:
            using Communicator = Traits<Mesh>::Communicator;

            class Impl;

            Mesh();
            ~Mesh();

            const Communicator &comm() const;
            Communicator &comm();

            bool read(const Path &path);
            void read(Input &in) override;
            void describe(std::ostream &os) const override;

        private:
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace sfem
}  // namespace utopia

#endif  // UTOPIA_SFEM_MESH_HPP

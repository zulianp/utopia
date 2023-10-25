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
            using Traits = utopia::Traits<Mesh>;
            using Communicator = Traits::Communicator;
            using Vector = Traits::Vector;
            using SizeType = Traits::SizeType;

            class Impl;

            Mesh();
            ~Mesh();

            const Communicator &comm() const;
            Communicator &comm();

            bool read(const Path &path);
            bool write(const Path &path);
            void read(Input &in) override;
            void write_nodal_field(const Path &path, const Vector &field);
            void describe(std::ostream &os) const override;

            void create_vector_nodal(Vector &out, int components = 1) const;

            SizeType n_local_nodes() const;

            ArrayView<const SizeType> node_mapping() const;

            void *raw_type() const;

        private:
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace sfem
}  // namespace utopia

#endif  // UTOPIA_SFEM_MESH_HPP

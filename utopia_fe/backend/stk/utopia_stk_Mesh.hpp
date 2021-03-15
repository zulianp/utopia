#ifndef UTOPIA_STK_MESH_HPP
#define UTOPIA_STK_MESH_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
// #include "utopia_Mesh.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_stk_ForwardDeclarations.hpp"

#include <memory>

namespace utopia {

    template <>
    class Traits<utopia::stk::Mesh> : public Traits<UVector> {};

    namespace stk {

        class Mesh final : public Configurable, public Describable {
        public:
            using SizeType = Traits<Mesh>::SizeType;
            using Comm = Traits<Mesh>::Communicator;

            ~Mesh();
            Mesh(const Comm &comm = Comm::get_default());

            bool read(const Path &path);
            bool write(const Path &path);

            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            const Comm &comm() const;

            ::stk::mesh::BulkData &raw_type();
            const ::stk::mesh::BulkData &raw_type() const;
            void wrap(const std::shared_ptr<::stk::mesh::BulkData> &mesh);
            bool empty() const;

            void unit_cube(const SizeType &nx = 10, const SizeType &ny = 10, const SizeType &nz = 10);

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_MESH_HPP

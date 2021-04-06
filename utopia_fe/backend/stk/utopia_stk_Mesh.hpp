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
            using Vector = Traits<Mesh>::Vector;
            using Matrix = Traits<Mesh>::Matrix;

            using Comm = Traits<Mesh>::Communicator;

            ~Mesh();
            Mesh(const Comm &comm = Comm::get_default());

            bool read(const Path &path);
            bool write(const Path &path) const;

            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            const Comm &comm() const;

            ::stk::mesh::BulkData &bulk_data() const;
            ::stk::mesh::MetaData &meta_data() const;

            void wrap(const std::shared_ptr<::stk::mesh::MetaData> &meta_data,
                      const std::shared_ptr<::stk::mesh::BulkData> &bulk_data);

            bool empty() const;
            int spatial_dimension() const;

            SizeType n_elements() const;
            SizeType n_nodes() const;

            SizeType n_local_elements() const;
            SizeType n_local_nodes() const;

            void displace(const Vector &displacement);

            void unit_cube(const SizeType &nx, const SizeType &ny, const SizeType &nz);

        private:
            class Impl;
            class IO;
            std::unique_ptr<Impl> impl_;

            friend class IO;
            friend class Impl;
        };

    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_MESH_HPP

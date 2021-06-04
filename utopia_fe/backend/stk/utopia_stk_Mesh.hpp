#ifndef UTOPIA_STK_MESH_HPP
#define UTOPIA_STK_MESH_HPP

#include "utopia_AABB.hpp"
#include "utopia_Communicator.hpp"
#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_stk_ForwardDeclarations.hpp"

#include "utopia_stk_SideSet.hpp"

#include <memory>

namespace utopia {

    template <>
    class Traits<utopia::stk::Mesh> : public Traits<UVector> {
    public:
        using Super = Traits<UVector>;
        using SideSet = utopia::stk::SideSet;
        using AABB = utopia::AABB<std::vector<Super::Scalar>>;
    };

    namespace stk {

        class Mesh final : public Configurable, public Describable {
        public:
            using SizeType = Traits<Mesh>::SizeType;
            using Scalar = Traits<Mesh>::Scalar;
            using Vector = Traits<Mesh>::Vector;
            using Matrix = Traits<Mesh>::Matrix;
            using AABB = Traits<Mesh>::AABB;

            using Comm = Traits<Mesh>::Communicator;

            ~Mesh();
            Mesh(const Comm &comm = Comm::get_default());

            bool read(const Path &path);
            bool write(const Path &path);

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
            void scale(const Scalar &scale_factor);

            void unit_cube(const SizeType &nx, const SizeType &ny, const SizeType &nz);

            void box(const AABB &box,
                     const SizeType &nx,
                     const SizeType &ny,
                     const SizeType &nz,
                     const std::string &elem_type = "HEX8");

            void bounding_box(AABB &output) const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;

            friend class MeshIO;
            friend class Impl;

            void init();
        };

    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_MESH_HPP

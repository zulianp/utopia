#ifndef UTOPIA_MARS_MESH_HPP
#define UTOPIA_MARS_MESH_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_mars_ForwardDeclarations.hpp"

// FIXME
#include "mars_context.hpp"

#include <memory>

namespace utopia {

    template <>
    class Traits<utopia::mars::Mesh> : public Traits<TpetraVector> {
    public:
        // using SideSet = utopia::mars::SideSet;
    };

    namespace mars {

        class Mesh final : public Configurable, public Describable {
        public:
            using SizeType = Traits<Mesh>::SizeType;
            using Scalar = Traits<Mesh>::Scalar;
            using Vector = Traits<Mesh>::Vector;
            using Matrix = Traits<Mesh>::Matrix;

            using Comm = Traits<Mesh>::Communicator;

            ~Mesh();
            Mesh(const Comm &comm = Comm::get_default());

            // bool read(const Path &path);
            // bool write(const Path &path);

            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            const Comm &comm() const;

            bool empty() const;
            int spatial_dimension() const;
            int manifold_dimension() const;

            SizeType n_elements() const;
            SizeType n_nodes() const;

            SizeType n_local_elements() const;
            SizeType n_local_nodes() const;

            /// Users should not try to use this
            template <class RawType>
            std::shared_ptr<RawType> raw_type() const;

            ::mars::context &raw_type_context();

            // void displace(const Vector &displacement);
            // void scale(const Scalar &scale_factor);

            void unit_cube(const SizeType &nx, const SizeType &ny, const SizeType &nz);

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;

            void init();
        };

    }  // namespace mars

}  // namespace utopia

#endif  // UTOPIA_MARS_MESH_HPP

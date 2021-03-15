#ifndef UTOPIA_MOONOLITH_MESH_HPP
#define UTOPIA_MOONOLITH_MESH_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_moonolith_ForwardDeclarations.hpp"

#include <memory>

namespace utopia {

    template <>
    class Traits<utopia::moonolith::Mesh> : public Traits<UVector> {};

    namespace moonolith {

        class Mesh final : public Configurable, public Describable {
        public:
            using SizeType = Traits<Mesh>::SizeType;
            using Scalar = Traits<Mesh>::Scalar;
            using Comm = Traits<Mesh>::Communicator;

            ~Mesh();
            Mesh(const Comm &comm = Comm::get_default());

            bool read(const Path &path);
            bool write(const Path &path) const;

            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            const Comm &comm() const;

            template <int Dim>
            ::moonolith::Mesh<Scalar, Dim> &raw_type();

            template <int Dim>
            const ::moonolith::Mesh<Scalar, Dim> &raw_type() const;

            template <int Dim>
            void wrap(const std::shared_ptr<::moonolith::Mesh<Scalar, Dim>> &mesh);

            bool empty() const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace moonolith
}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_MESH_HPP
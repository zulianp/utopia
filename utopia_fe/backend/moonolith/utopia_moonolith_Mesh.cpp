#include "utopia_moonolith_Mesh.hpp"
#include "par_moonolith.hpp"

namespace utopia {

    namespace moonolith {
        class Mesh::Impl {};

        Mesh::~Mesh() {}

        Mesh::Mesh(const Comm &comm) {}

        bool Mesh::read(const Path &path) {}

        bool Mesh::write(const Path &path) const {}

        void Mesh::read(Input &in) {}

        void Mesh::describe(std::ostream &os) const {}

        const Mesh::Comm &Mesh::comm() const {}

        template <int Dim>
        ::moonolith::Mesh<Mesh::Scalar, Dim> &Mesh::raw_type() {}

        template <int Dim>
        const ::moonolith::Mesh<Mesh::Scalar, Dim> &Mesh::raw_type() const {}

        template <int Dim>
        void Mesh::wrap(const std::shared_ptr<::moonolith::Mesh<Scalar, Dim>> &mesh) {}

        bool Mesh::empty() const {}

    }  // namespace moonolith
}  // namespace utopia
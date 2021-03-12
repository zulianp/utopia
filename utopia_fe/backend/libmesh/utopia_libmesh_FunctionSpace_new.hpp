#ifndef UTOPIA_LIBMESH_FUNCTION_SPACE_NEW_HPP
#define UTOPIA_LIBMESH_FUNCTION_SPACE_NEW_HPP

#include "utopia_libmesh_Mesh.hpp"

namespace utopia {
    template <>
    class Traits<utopia::libmesh::FunctionSpace> : public Traits<utopia::libmesh::Mesh> {};

    namespace libmesh {

        class FunctionSpace : public Configurable, public Describable {
        public:
            FunctionSpace(const Communicator &comm = Traits<FunctionSpace>::Communicator::get_default());
            FunctionSpace(const std::shared_ptr<Mesh> &mesh);
            ~FunctionSpace();

            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            std::shared_ptr<Mesh> mesh_ptr() const;
            const Mesh &mesh() const;
            Mesh &mesh();

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;

            class Sys;
            class Var;
            class BC;
        };
    }  // namespace libmesh

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_FUNCTION_SPACE_NEW_HPP

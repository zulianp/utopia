#ifndef UTOPIA_LIBMESH_FUNCTION_SPACE_HPP
#define UTOPIA_LIBMESH_FUNCTION_SPACE_HPP

// utopia
#include "utopia_FunctionSpace.hpp"

// utopia_libmesh
#include "utopia_libmesh_Mesh.hpp"

namespace utopia {

    template <>
    class FunctionSpace<LMMesh> final : public IFunctionSpace {
    public:
        FunctionSpace(const Communicator &comm);
        ~FunctionSpace();

        void read(Input &is) override;
        void describe(std::ostream &os = std::cout) const override;

    private:
        std::shared_ptr<LMMesh> mesh_;
        std::shared_ptr<libMesh::EquationSystems> equation_systems_;

        class Impl;
        std::unique_ptr<Impl> impl_;
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_FUNCTION_SPACE_HPP

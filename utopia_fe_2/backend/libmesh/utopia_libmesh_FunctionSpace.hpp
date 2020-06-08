#ifndef UTOPIA_LIBMESH_FUNCTION_SPACE_HPP
#define UTOPIA_LIBMESH_FUNCTION_SPACE_HPP

// utopia
#include "utopia_FunctionSpace.hpp"

// utopia_libmesh
#include "utopia_libmesh_Mesh.hpp"

namespace utopia {

    template <>
    class FunctionSpace<LMMesh> final : public IFunctionSpace<LMMesh> {
    public:
        using Traits = utopia::Traits<LMMesh>;
        using Vector = typename Traits::Vector;
        using Matrix = typename Traits::Matrix;
        using Communicator = typename Traits::Communicator;

        FunctionSpace(const Communicator &comm);
        ~FunctionSpace();

        void read(Input &is) override;
        void describe(std::ostream &os = std::cout) const override;
        bool write(const Path &path, const Vector &x) const override;
        void create_vector(Vector &x) const override;

    private:
        std::shared_ptr<LMMesh> mesh_;
        std::shared_ptr<libMesh::EquationSystems> equation_systems_;

        class Impl;
        std::unique_ptr<Impl> impl_;

        libMesh::System &system();
        const libMesh::System &system() const;
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_FUNCTION_SPACE_HPP

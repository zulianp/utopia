#ifndef UTOPIA_LIBMESH_FUNCTION_SPACE_HPP
#define UTOPIA_LIBMESH_FUNCTION_SPACE_HPP

// utopia
#include "utopia_FunctionSpace.hpp"
#include "utopia_blas_Matrix.hpp"
#include "utopia_blas_Vector.hpp"

// utopia_libmesh
#include "utopia_libmesh_Mesh.hpp"

#include "libmesh/fe_base.h"

namespace utopia {

    template <>
    class FunctionSpace<LMMesh> final : public IFunctionSpace<LMMesh> {
    public:
        using Traits = utopia::Traits<LMMesh>;
        using SizeType = Traits::SizeType;
        using Scalar = Traits::Scalar;
        using Vector = Traits::Vector;
        using Matrix = Traits::Matrix;
        using Communicator = Traits::Communicator;
        using Elem = libMesh::Elem;
        using ElementMatrix = utopia::BlasMatrix<Scalar>;
        using ElementVector = utopia::BlasVector<Scalar>;

        FunctionSpace(const Communicator &comm);
        ~FunctionSpace();

        void read(Input &is) override;
        void describe(std::ostream &os = std::cout) const override;
        bool write(const Path &path, const Vector &x) const override;

        void create_vector(Vector &x) const override;
        void create_matrix(Matrix &A) const override;

        void apply_constraints(Vector &x) const override;
        void apply_constraints(Matrix &A, Vector &b) const override;
        void apply_zero_constraints(Vector &x) const override;

        void add_matrix(const Elem &e, const ElementMatrix &el_mat, Matrix &mat) const;
        void add_vector(const Elem &e, const ElementVector &el_vec, Vector &vec) const;

        inline LMMesh &mesh() { return *mesh_; }
        inline const LMMesh &mesh() const { return *mesh_; }

        int dim() const;
        int order(const int var_num) const;
        std::unique_ptr<libMesh::FEBase> make_fe(const int var_num) const;

        libMesh::DofMap &dof_map();
        const libMesh::DofMap &dof_map() const;

    private:
        std::shared_ptr<LMMesh> mesh_;
        std::shared_ptr<libMesh::EquationSystems> equation_systems_;

        class Impl;
        std::unique_ptr<Impl> impl_;

        libMesh::System &system();
        const libMesh::System &system() const;
    };

    using LMFunctionSpace = utopia::FunctionSpace<LMMesh>;

    template <>
    class Traits<FunctionSpace<LMMesh>> : public Traits<LMMesh> {
    public:
        using ElementMatrix = utopia::BlasMatrix<Scalar>;
        using ElementVector = utopia::BlasVector<Scalar>;
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_FUNCTION_SPACE_HPP

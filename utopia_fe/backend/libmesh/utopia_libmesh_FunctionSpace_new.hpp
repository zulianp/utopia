#ifndef UTOPIA_LIBMESH_FUNCTION_SPACE_NEW_HPP
#define UTOPIA_LIBMESH_FUNCTION_SPACE_NEW_HPP

#include "utopia_Field.hpp"
#include "utopia_FunctionSpaceBase.hpp"

#include "utopia_libmesh_Mesh.hpp"
#include "utopia_libmesh_SideSet.hpp"

namespace utopia {
    template <>
    class Traits<utopia::libmesh::FunctionSpace> : public Traits<utopia::libmesh::Mesh> {
    public:
        using Mesh = utopia::libmesh::Mesh;
    };

    template <>
    class Traits<utopia::libmesh::FunctionSubspace> : public Traits<utopia::libmesh::FunctionSpace> {};

    namespace libmesh {

        class FunctionSpaceWrapper;

        class FunctionSubspace : public Traits<FunctionSubspace> {
        public:
            FunctionSubspace();
            ~FunctionSubspace();

            friend class utopia::libmesh::FunctionSpace;

            SizeType n_dofs() const;
            SizeType n_local_dofs() const;
            SizeType n_subspaces() const;
            SizeType subspace_id() const;

            // access main function space subspaces (main system)
            FunctionSubspace subspace(const SizeType i, const SizeType n_vars = 1);
            inline FunctionSubspace operator[](const SizeType i) { return subspace(i); }

        private:
            class Impl;
            std::shared_ptr<Impl> impl_;
        };

        class FunctionSpace : public FunctionSpaceBase<libmesh::Mesh> {
        public:
            using Vector = Traits<FunctionSpace>::Vector;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using SizeType = Traits<FunctionSpace>::SizeType;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using IndexArray = Traits<FunctionSpace>::IndexArray;
            using Comm = Traits<FunctionSpace>::Communicator;

            FunctionSpace(const Comm &comm = Comm::get_default());
            FunctionSpace(const std::shared_ptr<Mesh> &mesh);
            void init(const std::shared_ptr<Mesh> &mesh) override;
            void update(const SimulationTime<Scalar> &) override;

            ~FunctionSpace();

            bool write(const Path &path, const Vector &x) override;
            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            std::shared_ptr<Mesh> mesh_ptr() const override;
            const Mesh &mesh() const override;
            Mesh &mesh() override;

            inline const Comm &comm() const override { return mesh().comm(); }

            libMesh::EquationSystems &raw_type();
            const libMesh::EquationSystems &raw_type() const;
            const libMesh::DofMap &raw_type_dof_map() const;
            libMesh::DofMap &raw_type_dof_map();
            const libMesh::System &raw_type_system() const;

            // The binary is private, so do not try to use it
            std::shared_ptr<FunctionSpaceWrapper> wrapper();

            SizeType n_dofs() const override;
            SizeType n_local_dofs() const override;
            SizeType n_subspaces() const;
            SizeType system_id() const;

            int n_var() const override;

            // access main function space subspaces (main system)
            FunctionSubspace subspace(const SizeType i, const SizeType n_vars = 1);
            inline FunctionSubspace operator[](const SizeType i) { return subspace(i); }

            // For aux systems
            SizeType n_systems() const;
            // access other spaces (auxiliry systems)
            FunctionSubspace auxiliary_space(const SizeType i);

            void node_eval(std::function<void(const SizeType idx, const Scalar *)> fun);

            void create_matrix(Matrix &mat) const override;
            void create_vector(Vector &vec) const override;
            void create_field(Field<FunctionSpace> &field);

            void apply_constraints(Matrix &mat, Vector &vec) const override;
            void apply_constraints(Vector &vec) const override;
            void apply_constraints(Matrix &mat, const Scalar diag_value = 1) const override;

            void create_boundary_node_list(IndexArray &indices) const;

            void apply_zero_constraints(Vector &vec) const override;
            void add_dirichlet_boundary_condition(const std::string &boundary_name,
                                                  const Scalar &value,
                                                  const int variable = 0) override;

            void add_dirichlet_boundary_condition(const int boundary_id, const Scalar &value, const int variable = 0);

            bool empty() const override;

            bool read(const Path &path,
                      const std::vector<std::string> &var_names,
                      Vector &val,
                      const int time_step = 1);

            bool read_with_state(Input &in, Field<FunctionSpace> &val);

            void displace(const Vector &displacement) override;

            const std::string &name() const override;

            void copy_meta_info_from(const FunctionSpace &other);

            void initialize() override;

        private:
            using Impl = utopia::libmesh::FunctionSpaceWrapper;
            std::shared_ptr<Impl> impl_;

            class Sys;
            class Var;
            class BC;
        };

    }  // namespace libmesh

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_FUNCTION_SPACE_NEW_HPP

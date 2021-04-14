#ifndef UTOPIA_LIBMESH_FUNCTION_SPACE_NEW_HPP
#define UTOPIA_LIBMESH_FUNCTION_SPACE_NEW_HPP

#include "utopia_Field.hpp"
#include "utopia_libmesh_Mesh.hpp"

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

        class FunctionSpace : public Configurable, public Describable, public Traits<FunctionSpace> {
        public:
            using Vector = Traits<FunctionSpace>::Vector;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Comm = Traits<FunctionSpace>::Communicator;

            FunctionSpace(const Comm &comm = Comm::get_default());
            FunctionSpace(const std::shared_ptr<Mesh> &mesh);
            ~FunctionSpace();

            bool write(const Path &path, const Vector &x);
            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            std::shared_ptr<Mesh> mesh_ptr() const;
            const Mesh &mesh() const;
            Mesh &mesh();

            inline const Comm &comm() const { return mesh().comm(); }

            libMesh::EquationSystems &raw_type();
            const libMesh::EquationSystems &raw_type() const;
            const libMesh::DofMap &raw_type_dof_map() const;
            const libMesh::System &raw_type_system() const;

            // The binary is private, so do not try to use it
            std::shared_ptr<FunctionSpaceWrapper> wrapper();

            SizeType n_dofs() const;
            SizeType n_local_dofs() const;
            SizeType n_subspaces() const;
            SizeType system_id() const;

            // access main function space subspaces (main system)
            FunctionSubspace subspace(const SizeType i, const SizeType n_vars = 1);
            inline FunctionSubspace operator[](const SizeType i) { return subspace(i); }

            // For aux systems
            SizeType n_systems() const;
            // access other spaces (auxiliry systems)
            FunctionSubspace auxiliary_space(const SizeType i);

            void create_matrix(Matrix &mat) const;
            void create_vector(Vector &vec) const;
            void create_field(Field<FunctionSpace> &field);

            void apply_constraints(Matrix &mat, Vector &vec) const;
            void apply_constraints(Vector &vec) const;
            void apply_constraints(Matrix &mat) const;

            void apply_zero_constraints(Vector &vec) const;

            bool empty() const;

            bool read(const Path &path,
                      const std::vector<std::string> &var_names,
                      Vector &val,
                      const int time_step = 1);

            bool read_with_state(Input &in, Field<FunctionSpace> &val);

            void displace(const Vector &displacement);

            const std::string &name() const;

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

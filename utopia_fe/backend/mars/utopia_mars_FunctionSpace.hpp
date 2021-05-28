#ifndef UTOPIA_MARS_FUNCTION_SPACE_HPP
#define UTOPIA_MARS_FUNCTION_SPACE_HPP

#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_mars_Mesh.hpp"

namespace utopia {
    template <>
    class Traits<utopia::mars::FunctionSpace> : public Traits<utopia::mars::Mesh> {
    public:
        using Mesh = utopia::mars::Mesh;
    };

    namespace mars {

        class FunctionSpace : public Configurable, public Describable, public Traits<FunctionSpace> {
        public:
            using Vector = Traits<FunctionSpace>::Vector;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;
            using LocalSizeType = Traits<Matrix>::LocalSizeType;
            using IndexSet = Traits<FunctionSpace>::IndexSet;
            using Comm = Traits<FunctionSpace>::Communicator;

            FunctionSpace(const Comm &comm = Comm::get_default());
            FunctionSpace(const std::shared_ptr<Mesh> &mesh);
            ~FunctionSpace();

            void init(const std::shared_ptr<Mesh> &mesh);

            bool write(const Path &path, const Vector &x);
            void read(Input &in) override;
            // bool read_with_state(Input &in, Field<FunctionSpace> &val);
            void describe(std::ostream &os = std::cout) const override;

            std::shared_ptr<Mesh> mesh_ptr() const;
            const Mesh &mesh() const;
            Mesh &mesh();

            inline const Comm &comm() const { return mesh().comm(); }

            SizeType n_dofs() const;
            SizeType n_local_dofs() const;
            int n_var() const;
            void set_n_var(const int n_var);

            void create_vector(Vector &v) const;
            void create_local_vector(Vector &v) const;
            void create_matrix(Matrix &m) const;

            // void create_field(Field<FunctionSpace> &field);
            // void create_nodal_vector_field(const int vector_size, Field<FunctionSpace> &field);

            void apply_constraints(Matrix &m, const Scalar diag_value = 1.0);
            void apply_constraints(Vector &v);
            void apply_constraints(Matrix &m, Vector &v);
            void apply_zero_constraints(Vector &vec) const;

            void add_dirichlet_boundary_condition(const std::string &name,
                                                  const Scalar &value,
                                                  const int component = 0);

            bool empty() const;

            // void displace(const Vector &displacement);
            void global_to_local(const Vector &global, Vector &local) const;
            void local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const;

            // const DofMap &dof_map() const;
            // DofMap &dof_map();

            const std::string &name() const;

            // SizeType n_variables() const;
            // const std::string &variable_name(const SizeType var_num) const;
            // int variable_size(const SizeType var_num) const;

            // void nodal_field_to_local_vector(Vector &v);
            // void local_vector_to_nodal_field(const Vector &v);

            // void nodal_field_to_global_vector(Vector &v);
            // void global_vector_to_nodal_field(const Vector &v);

            // template <typename FieldType>
            // void declare_new_nodal_field(const std::string &name, const int n_comp);

            // void create_vector(const std::string &field_name, Vector &v) const;
            // void global_vector_to_nodal_field(const std::string &field_name, const Vector &v);

            // void backend_set_nodal_field(const Field<FunctionSpace> &field);

        private:
            class Impl;
            // class Var;

            std::shared_ptr<Impl> impl_;

            // void register_output_variables(MeshIO &io);
            // void read_meta(Input &in);
            // void register_variables();
            // friend class SpaceIO;
        };

    }  // namespace mars

}  // namespace utopia

#endif  // UTOPIA_MARS_FUNCTION_SPACE_HPP

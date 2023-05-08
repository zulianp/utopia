#ifndef UTOPIA_STK_FUNCTION_SPACE_HPP
#define UTOPIA_STK_FUNCTION_SPACE_HPP

#include "utopia_DirichletBoundary.hpp"
#include "utopia_FEVar.hpp"
#include "utopia_FunctionSpaceBase.hpp"
#include "utopia_fe_Core.hpp"

#include "utopia_Field.hpp"
#include "utopia_stk_Mesh.hpp"

namespace utopia {
    template <>
    class Traits<utopia::stk::FunctionSpace> : public Traits<utopia::stk::Mesh> {
    public:
        using Mesh = utopia::stk::Mesh;
        using Environment = utopia::Environment<utopia::stk::FunctionSpace>;
    };

    namespace stk {

        class FunctionSpace : public FunctionSpaceBase<stk::Mesh>, public Traits<FunctionSpace> {
        public:
            using Vector = Traits<FunctionSpace>::Vector;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;
            using IndexSet = Traits<FunctionSpace>::IndexSet;
            using IndexArray = Traits<FunctionSpace>::IndexArray;
            using Comm = Traits<FunctionSpace>::Communicator;
            using DirichletBoundary = utopia::DirichletBoundary<Traits<FunctionSpace>>;

            FunctionSpace(const Comm &comm = Comm::get_default());
            FunctionSpace(const std::shared_ptr<Mesh> &mesh);
            ~FunctionSpace();

            void init(const std::shared_ptr<Mesh> &mesh) override;
            void update(const SimulationTime<Scalar> &) override;

            bool write(const Path &path, const Vector &x) override;
            void read(Input &in) override;
            bool read_with_state(Input &in, Field<FunctionSpace> &val);
            void describe(std::ostream &os) const override;

            std::shared_ptr<Mesh> mesh_ptr() const override;
            const Mesh &mesh() const override;
            Mesh &mesh() override;

            inline const Comm &comm() const override { return mesh().comm(); }

            SizeType n_dofs() const override;
            SizeType n_local_dofs() const override;
            int n_var() const override;
            void set_n_var(const int n_var);

            void create_vector(Vector &v) const override;

            void create_local_vector(Vector &v) const;
            void create_matrix(Matrix &m) const override;
            void create_field(Field<FunctionSpace> &field);
            void create_nodal_vector_field(const int vector_size, Field<FunctionSpace> &field);

            void apply_constraints(Matrix &m, const Scalar diag_value = 1.0) const override;
            void apply_constraints(Vector &v) const override;
            void apply_constraints(Matrix &m, Vector &v) const override;
            void apply_zero_constraints(Vector &vec) const override;
            void copy_at_constrained_nodes(const Vector &, Vector &) const override;

            void overwrite_parts(const std::vector<std::string> &parts,
                                 const std::vector<int> &components,
                                 const Vector &source,
                                 Vector &destination) const;

            void set_overwrite_vector(const Vector &v);

            void add_dirichlet_boundary_condition(const std::string &name,
                                                  const Scalar &value,
                                                  const int component = 0) override;

            bool empty() const override;

            // void displacement_field_from_transform(const std::vector<Scalar> &scale_factors,
            //                                        Field<FunctionSpace> &displacement);

            void displace(const Vector &displacement) override;
            void global_to_local(const Vector &global, Vector &local) const;
            void local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const;

            const DofMap &dof_map() const;
            DofMap &dof_map();

            const std::string &name() const override;

            SizeType n_variables() const;
            const std::string &variable_name(const SizeType var_num) const;
            int variable_size(const SizeType var_num) const;

            int add_variable(const FEVar &var);

            void nodal_field_to_local_vector(Vector &v);
            void local_vector_to_nodal_field(const Vector &v);

            void nodal_field_to_global_vector(Vector &v);
            void global_vector_to_nodal_field(const Vector &v);

            template <typename FieldType>
            void declare_new_nodal_field(const std::string &name, const int n_comp);

            void create_vector(const std::string &field_name, Vector &v) const;
            void global_vector_to_nodal_field(const std::string &field_name, const Vector &v);

            void backend_set_nodal_field(const Field<FunctionSpace> &field);
            void backend_set_elemental_field(const Field<FunctionSpace> &field);

            void copy_meta_info_from(const FunctionSpace &other);
            void initialize() override { initialize(true); }
            void initialize(const bool valid_local_id_mode);

            void node_eval(std::function<void(const SizeType idx, const Scalar *)> fun);
            void node_eval(const std::string &part_name, std::function<void(const SizeType idx, const Scalar *)> fun);

            const DirichletBoundary &dirichlet_boundary() const;

            void set_print_map(const bool val);

            void create_boundary_node_list(IndexArray &node_list) const;

            void create_node_to_element_matrix(Matrix &matrix) const;

        private:
            class Impl;
            // class Var;

            std::shared_ptr<Impl> impl_;

            void register_output_variables(MeshIO &io);
            void read_meta(Input &in);
            void register_variables();
            friend class SpaceIO;
        };

    }  // namespace stk

}  // namespace utopia

#endif  // UTOPIA_STK_FUNCTION_SPACE_HPP

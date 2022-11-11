#ifndef UTOPIA_MARS_FUNCTION_SPACE_HPP
#define UTOPIA_MARS_FUNCTION_SPACE_HPP

#include <memory>
#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"

#include "utopia_DirichletBoundary.hpp"

#include "utopia_mars_ForwardDeclarations.hpp"
#include "utopia_mars_Mesh.hpp"

namespace utopia {
    template <>
    class Traits<utopia::mars::FunctionSpace> : public Traits<utopia::mars::Mesh> {
    public:
        using Mesh = utopia::mars::Mesh;
        using Environment = utopia::Environment<utopia::mars::FunctionSpace>;
    };

    namespace mars {

        class IFEHandler {
        public:
            using Vector = Traits<FunctionSpace>::Vector;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;
            using LocalSizeType = Traits<Matrix>::LocalSizeType;
            using IndexSet = Traits<FunctionSpace>::IndexSet;
            using Comm = Traits<FunctionSpace>::Communicator;

            using MarsCrsMatrix = Matrix::CrsMatrixType::local_matrix_type;

            virtual ~IFEHandler() = default;

            virtual MarsCrsMatrix new_crs_matrix() = 0;

            virtual void describe() const = 0;
            virtual void matrix_apply_constraints(Matrix &m,
                                                  const Scalar diag_value,
                                                  const std::string side,
                                                  const int component) = 0;

            virtual void vector_apply_constraints(Vector &v,
                                                  const Scalar value,
                                                  const std::string side,
                                                  const int component) = 0;

            virtual void apply_zero_constraints(Vector &vec, const std::string side, const int component) = 0;

            virtual void copy_at_constrained_nodes(const Vector &in,
                                                   Vector &out,
                                                   const std::string side,
                                                   const int component) = 0;

            // virtual void system_apply_constraints(Matrix &m, Vector &v) = 0;

            virtual void ensure_sparsity_pattern() = 0;

            virtual SizeType n_local_dofs() = 0;

            virtual SizeType n_dofs() = 0;

            virtual Factory &factory() = 0;
        };

        class FunctionSpace : public Configurable, public Describable, public Traits<FunctionSpace> {
        public:
            using Vector = Traits<FunctionSpace>::Vector;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;
            using LocalSizeType = Traits<Matrix>::LocalSizeType;
            using IndexSet = Traits<FunctionSpace>::IndexSet;
            using Comm = Traits<FunctionSpace>::Communicator;
            using DirichletBoundary = utopia::DirichletBoundary<Traits<FunctionSpace>>;

            FunctionSpace(const Comm &comm = Comm::get_default());
            FunctionSpace(const std::shared_ptr<Mesh> &mesh);
            ~FunctionSpace();

            void init(const std::shared_ptr<Mesh> &mesh);
            void update(const SimulationTime<Scalar> &);

            bool write(const Path &path, const Vector &x);
            void read(Input &in) override;
            // bool read_with_state(Input &in, Field<FunctionSpace> &val);
            void describe(std::ostream &os = std::cout) const override;

            std::shared_ptr<Mesh> mesh_ptr() const;
            const Mesh &mesh() const;
            Mesh &mesh();

            template <class RawType>
            std::shared_ptr<RawType> raw_type() const;

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

            void apply_constraints(Matrix &m, const Scalar diag_value = 1.0) const;
            void apply_constraints(Vector &v) const;
            void apply_constraints(Matrix &m, Vector &v) const;
            void apply_zero_constraints(Vector &vec) const;
            void copy_at_constrained_nodes(const Vector &in, Vector &out) const /*override*/;

            void apply_constraints_update(Vector &v) const { this->apply_constraints(v); }

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

            const Factory &factory() const;

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
            std::unique_ptr<Impl> impl_;
            std::shared_ptr<IFEHandler> handler() const;

            // void register_output_variables(MeshIO &io);
            // void read_meta(Input &in);
            // void register_variables();
            // friend class SpaceIO;
        };

    }  // namespace mars

}  // namespace utopia

#endif  // UTOPIA_MARS_FUNCTION_SPACE_HPP

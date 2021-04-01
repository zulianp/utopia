#ifndef UTOPIA_STK_FUNCTION_SPACE_HPP
#define UTOPIA_STK_FUNCTION_SPACE_HPP

#include "utopia_stk_Mesh.hpp"

namespace utopia {
    template <>
    class Traits<utopia::stk::FunctionSpace> : public Traits<utopia::stk::Mesh> {
    public:
        using Mesh = utopia::stk::Mesh;
    };

    namespace stk {

        class FunctionSpace : public Configurable, public Describable, public Traits<FunctionSpace> {
        public:
            using Vector = Traits<FunctionSpace>::Vector;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;
            using IndexSet = Traits<FunctionSpace>::IndexSet;
            using Comm = Traits<FunctionSpace>::Communicator;

            FunctionSpace(const Comm &comm = Comm::get_default());
            FunctionSpace(const std::shared_ptr<Mesh> &mesh);
            ~FunctionSpace();

            void init(const std::shared_ptr<Mesh> &mesh);

            bool write(const Path &path, const Vector &x);
            void read(Input &in) override;
            void read_with_state(Input &in, Vector &val);
            void describe(std::ostream &os) const override;

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

            void apply_constraints(Matrix &m);
            void apply_constraints(Vector &v);
            void apply_constraints(Matrix &m, Vector &v);

            void add_dirichlet_boundary_condition(const std::string &name,
                                                  const Scalar &value,
                                                  const int component = 0);

            bool empty() const;

            void displace(const Vector &displacement);
            void global_to_local(const Vector &global, Vector &local) const;

            const DofMap &dof_map() const;
            DofMap &dof_map();

        private:
            class Impl;
            class Var;

            std::shared_ptr<Impl> impl_;
        };

    }  // namespace stk

}  // namespace utopia

#endif  // UTOPIA_STK_FUNCTION_SPACE_HPP

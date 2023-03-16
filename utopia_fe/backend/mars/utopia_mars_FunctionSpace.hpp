#ifndef UTOPIA_MARS_FUNCTION_SPACE_HPP
#define UTOPIA_MARS_FUNCTION_SPACE_HPP

#include <memory>
#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"

#include "utopia_DirichletBoundary.hpp"

#include "utopia_kokkos_Discretization.hpp"
#include "utopia_kokkos_UniformFE.hpp"
#include "utopia_mars_ForwardDeclarations.hpp"
#include "utopia_mars_Mesh.hpp"

namespace utopia {
    template <>
    class Traits<utopia::mars::FunctionSpace> : public Traits<utopia::mars::Mesh> {
    public:
        using Super_ = Traits<utopia::mars::Mesh>;
        using Mesh = utopia::mars::Mesh;
        using Environment = utopia::Environment<utopia::mars::FunctionSpace>;
        using FE = utopia::kokkos::UniformFE<typename Super_::Scalar>;
        using Discretization = utopia::Discretization<utopia::mars::FunctionSpace, FE>;
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

            using FE = Traits<FunctionSpace>::FE;
            using KokkosDiscretization = utopia::kokkos::Discretization<FunctionSpace, FE>;
            using Part = KokkosDiscretization::Part;

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

            virtual void ensure_sparsity_pattern() = 0;

            virtual SizeType n_local_dofs() = 0;

            virtual SizeType n_dofs() = 0;

            virtual Factory &factory() = 0;

            ////////////////////////////////////////////////////////////////////////////////////

            virtual void create(std::vector<KokkosDiscretization::FE_ptr> &fe,
                                int order,
                                const KokkosDiscretization::Part &part = KokkosDiscretization::all()) = 0;

            virtual void create_on_boundary(std::vector<KokkosDiscretization::FE_ptr> &fe,
                                            int order,
                                            const KokkosDiscretization::Part &part = KokkosDiscretization::all()) = 0;

            ////////////////////////////////////////////////////////////////////////////////////

            virtual void convert_field(const Field<FunctionSpace> &in,
                                       std::vector<std::shared_ptr<KokkosDiscretization::FEField>> &out,
                                       const KokkosDiscretization::Part &part = KokkosDiscretization::all()) = 0;

            virtual void convert_field(const std::vector<std::shared_ptr<KokkosDiscretization::FEField>> &in,
                                       Field<FunctionSpace> &out,
                                       const KokkosDiscretization::Part &part = KokkosDiscretization::all()) = 0;

            ////////////////////////////////////////////////////////////////////////////////////

            virtual void global_to_local(const Vector &vector,
                                         std::vector<KokkosDiscretization::VectorAccumulator> &element_vectors,
                                         const KokkosDiscretization::Part &part = KokkosDiscretization::all(),
                                         const int comp = 0) = 0;

            // Local to global

            virtual void local_to_global(const std::vector<KokkosDiscretization::MatrixAccumulator> &acc,
                                         AssemblyMode mode,
                                         Matrix &mat,
                                         const KokkosDiscretization::Part &part = KokkosDiscretization::all()) = 0;

            virtual void local_to_global(const std::vector<KokkosDiscretization::VectorAccumulator> &acc,
                                         AssemblyMode mode,
                                         Vector &vec,
                                         const KokkosDiscretization::Part &part = KokkosDiscretization::all()) = 0;

            virtual void local_to_global(const Comm &comm,
                                         const std::vector<KokkosDiscretization::ScalarAccumulator> &acc,
                                         std::vector<Scalar> &scalars,
                                         const Part &part = KokkosDiscretization::all()) = 0;

            virtual void local_to_global_on_boundary(
                const std::vector<KokkosDiscretization::VectorAccumulator> &acc,
                AssemblyMode mode,
                Vector &vec,
                const KokkosDiscretization::Part &part = KokkosDiscretization::all()) = 0;
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

            using Discretization = Traits<FunctionSpace>::Discretization;

            using DirichletBoundary = utopia::DirichletBoundary<Traits<FunctionSpace>>;

            FunctionSpace(const Comm &comm = Comm::get_default());
            FunctionSpace(const std::shared_ptr<Mesh> &mesh);
            ~FunctionSpace();

            void init(const std::shared_ptr<Mesh> &mesh);
            void update(const SimulationTime<Scalar> &);

            bool write(const Path &path, const Vector &x);
            void read(Input &in) override;
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

            void global_to_local(const Vector &global, Vector &local) const;
            void local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const;

            const std::string &name() const;

            const Factory &factory() const;

            std::shared_ptr<IFEHandler> handler() const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace mars

}  // namespace utopia

#endif  // UTOPIA_MARS_FUNCTION_SPACE_HPP

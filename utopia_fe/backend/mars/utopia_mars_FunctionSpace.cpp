#include "utopia_mars_FunctionSpace.hpp"

#include "utopia_DirichletBoundary.hpp"
#include "utopia_FEVar.hpp"

#include "mars.hpp"
#include "mars_base.hpp"
#include "mars_context.hpp"
#include "mars_distributed_dof_management.hpp"
// #include "mars_distributed_finite_element.hpp"
// #include "mars_execution_context.hpp"
// #include "mars_imesh_kokkos.hpp"

namespace utopia {
    namespace mars {

        class FunctionSpace::Impl {
        public:
            // using DeviceType = typename Kokkos::Device<Kokkos::DefaultExecutionSpace, KokkosSpace>;
            using MarsCrsMatrix = Matrix::CrsMatrixType::local_matrix_type;
            using MapType = Matrix::MapType;

            using DeviceType = Matrix::CrsMatrixType::device_type;
            // using MarsCrsMatrix = typename KokkosSparse::CrsMatrix<Scalar, LocalSizeType, DeviceType, void,
            // SizeType>;

            template <class DMesh>
            void init(DMesh &mesh_impl) {
                static constexpr ::mars::Integer Degree = 1;
                using DofHandler = ::mars::DofHandler<DMesh, Degree>;
                using FEDofMap = ::mars::FEDofMap<DofHandler>;
                using SPattern =
                    ::mars::SparsityPattern<Scalar, LocalSizeType, SizeType, DofHandler, MarsCrsMatrix::size_type>;

                static_assert(std::is_same<SizeType, Matrix::CrsMatrixType::global_ordinal_type>::value, "Weird!");

                auto dof_handler_impl = std::make_shared<DofHandler>(&mesh_impl);  //, mesh->raw_type_context());
                dof_handler_impl->enumerate_dofs();

                dof_handler = dof_handler_impl;

                auto fe_dof_map_impl = std::make_shared<FEDofMap>(build_fe_dof_map(*dof_handler_impl));
                fe_dof_map = fe_dof_map_impl;

                SPattern sp(*dof_handler_impl);
                sp.build_pattern(*fe_dof_map_impl);

                new_crs_matrix = [sp]() -> MarsCrsMatrix { return sp.new_crs_matrix(); };

                describe = [dof_handler_impl, fe_dof_map_impl, sp]() {
                    dof_handler_impl->print_dofs();
                    sp.print_sparsity_pattern();

                    // ::mars::print_fe_dof_map(*dof_handler_impl, *fe_dof_map_impl);
                };

                matrix_apply_constraints = [sp](Matrix &m, const Scalar diag_value) {
                    // BC set constrained rows to zero, except diagonal where you set diag_value
                };

                vector_apply_constraints = [sp](Vector &v) {
                    // BC set values to constraint value (i.e., boundary value)
                };

                apply_zero_constraints = [sp](Vector &vec) {
                    // BC set values to constraint value to zero
                };

                system_apply_constraints = [this](Matrix &m, Vector &v) {
                    vector_apply_constraints(v);
                    matrix_apply_constraints(m, 1.0);
                };

                n_local_dofs = [dof_handler_impl]() { return dof_handler_impl->get_owned_dof_size(); };
                n_dofs = [dof_handler_impl]() { return dof_handler_impl->get_global_dof_size(); };
            }

            void read_meta(Input &in) {
                n_var = 0;

                if (!Options()
                         .add_option("name", name, "Unique name of the function space")
                         .add_option(
                             "n_var", n_var, "Number of variables per node (instead of specifiying all of them).")
                         .add_option("boundary_conditions", dirichlet_boundary, "Boundary conditions.")
                         .add_option("verbose", verbose, "Verbose output.")
                         .parse(in)) {
                    return;
                }

                in.get("variables", [this](Input &in) {
                    in.get_all([&](Input &in) {
                        FEVar v;
                        v.read(in);
                        variables.push_back(v);
                    });
                });

                if (variables.empty()) {
                    FEVar v;
                    // At least one variable
                    v.n_components = std::max(1, n_var);
                    variables.push_back(v);
                }

                int counted_vars = 0;

                for (auto &v : variables) {
                    counted_vars += v.n_components;
                }

                assert(n_var == 0 || n_var == counted_vars);

                n_var = counted_vars;
            }

            std::function<void()> describe;
            std::function<SizeType()> n_dofs;
            std::function<SizeType()> n_local_dofs;
            std::function<MarsCrsMatrix()> new_crs_matrix;

            std::function<void(Matrix &m, const Scalar diag_value)> matrix_apply_constraints;
            std::function<void(Vector &v)> vector_apply_constraints;
            std::function<void(Matrix &m, Vector &v)> system_apply_constraints;
            std::function<void(Vector &vec)> apply_zero_constraints;

            // std::function<crs_grap()> get_crs_graph;

            std::shared_ptr<::mars::IDofHandler> dof_handler;
            std::shared_ptr<::mars::IFEDofMap> fe_dof_map;
            std::shared_ptr<Mesh> mesh;

            std::string name;
            DirichletBoundary dirichlet_boundary;
            std::vector<FEVar> variables;

            // MarsCrsMatrix crs_matrix;

            bool verbose{false};
            int n_var{0};
        };

        FunctionSpace::FunctionSpace(const Comm &comm) : impl_(utopia::make_unique<Impl>()) {}

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>()) {
            init(mesh);
        }

        FunctionSpace::~FunctionSpace() {}

        void FunctionSpace::init(const std::shared_ptr<Mesh> &mesh) {
            impl_->mesh = mesh;
            switch (mesh->spatial_dimension()) {
                case 2: {
                    using DMesh = ::mars::DistributedMesh<::mars::ElementType::Quad4>;
                    auto mesh_impl = mesh->raw_type<DMesh>();
                    impl_->init(*mesh_impl);
                }
                case 3: {
                    using DMesh = ::mars::DistributedMesh<::mars::ElementType::Hex8>;
                    auto mesh_impl = mesh->raw_type<DMesh>();
                    impl_->init(*mesh_impl);
                }
            }
        }

        bool FunctionSpace::write(const Path &path, const Vector &x) {}

        void FunctionSpace::read(Input &in) {
            if (impl_->mesh) {
                impl_->mesh = std::make_shared<Mesh>();
                in.get("mesh", *impl_->mesh);
            }

            impl_->read_meta(in);
        }
        void FunctionSpace::describe(std::ostream &os) const {
            if (impl_->describe) {
                impl_->describe();
            } else {
                os << "FunctionSpace::describe: null\n";
            }
        }

        std::shared_ptr<FunctionSpace::Mesh> FunctionSpace::mesh_ptr() const { return impl_->mesh; }
        const FunctionSpace::Mesh &FunctionSpace::mesh() const { return *impl_->mesh; }
        FunctionSpace::Mesh &FunctionSpace::mesh() { return *impl_->mesh; }

        FunctionSpace::SizeType FunctionSpace::n_dofs() const { return impl_->n_dofs(); }
        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const { return impl_->n_local_dofs(); }
        int FunctionSpace::n_var() const { return 1; }
        void FunctionSpace::set_n_var(const int n_var) { assert(n_var == 1); }

        void FunctionSpace::create_vector(Vector &v) const { v.zeros(layout(comm(), n_local_dofs(), n_dofs())); }
        void FunctionSpace::create_local_vector(Vector &) const { Utopia::Abort("IMPLEMENT ME"); }

        void FunctionSpace::create_matrix(Matrix &m) const {
            using MapType = Matrix::MapType;

            SizeType n_global = this->n_dofs();
            SizeType n_local = this->n_local_dofs();

            const SizeType index_base = 0;
            ::Teuchos::RCP<MapType> row_map =
                rcp(new Matrix::MapType(n_global, n_local, index_base, this->comm().get()));
            // since MARS provides already global indices for the columns we can just use the identity map, replicated
            // on each process
            ::Teuchos::RCP<MapType> col_map =
                rcp(new Matrix::MapType(n_global, index_base, this->comm().get(), Tpetra::LocallyReplicated));

            // this constructor gives us a fill-complete matrix, so do not call fillComplete manually again
            Matrix::RCPCrsMatrixType mat =
                Teuchos::rcp(new Matrix::CrsMatrixType(impl_->new_crs_matrix(), row_map, col_map));
            m.wrap(mat, true);
        }

        void FunctionSpace::apply_constraints(Matrix &m, const Scalar diag_value) {}
        void FunctionSpace::apply_constraints(Vector &v) {}
        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) {}
        void FunctionSpace::apply_zero_constraints(Vector &vec) const {}

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name,
                                                             const Scalar &value,
                                                             const int component) {}

        bool FunctionSpace::empty() const { return !static_cast<bool>(impl_->new_crs_matrix); }

        void FunctionSpace::global_to_local(const Vector &global, Vector &local) const {}
        void FunctionSpace::local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const {}

        const std::string &FunctionSpace::name() const { return impl_->name; }

    }  // namespace mars

}  // namespace utopia

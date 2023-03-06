#include "utopia_mars_FunctionSpace.hpp"
#include <memory>

#include "utopia_DirichletBoundary.hpp"
#include "utopia_FEVar.hpp"
#include "utopia_mars_FEHandler.hpp"

#include "utopia_mars_MeshIO.hpp"

#include "utopia_mars_Discretization.hpp"

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
                using FEHandler = utopia::mars::FEHandler<DMesh>;

                auto handler_impl = std::make_shared<FEHandler>();
                assert(n_var != 0);
                handler_impl->init(mesh_impl, n_var);
                handler = handler_impl;

#ifdef MARS_WITH_WITH_IO
                write = [handler_impl, this](const Path &path, const Vector &x) -> bool {
                    MarsIOImpl<typename FEHandler::FEDofMap> w(handler_impl->get_fe_dof_map());

                    // auto x_kokkos = x.raw_type()->getLocalViewHost();
                    // return w.write_tpetra(path.to_string(), x_kokkos);

                    auto x_kokkos = x.raw_type()->getLocalViewHost();
                    w.add_field_tpetra("U", this->n_var, x_kokkos);
                    w.set_output_path(path);
                    return w.write();
                };
#else
                write = [](const Path &, const Vector &) -> bool { return false; };
#endif  // UTOPIA_ENABLE_VTK
            }

            template <class RawType>
            std::shared_ptr<RawType> raw_type() const {
                return std::dynamic_pointer_cast<RawType>(handler);
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

            std::shared_ptr<IFEHandler> handler;
            std::shared_ptr<Mesh> mesh;

            std::string name;
            DirichletBoundary dirichlet_boundary;
            std::vector<FEVar> variables;

            std::function<bool(const Path &, const Vector &)> write;

            bool verbose{false};
            int n_var{1};
        };

        std::shared_ptr<IFEHandler> FunctionSpace::handler() const { return impl_->handler; }

        FunctionSpace::FunctionSpace(const Comm &comm) : impl_(utopia::make_unique<Impl>()) {}

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>()) {
            init(mesh);
        }

        FunctionSpace::~FunctionSpace() {}

        void FunctionSpace::init(const std::shared_ptr<Mesh> &mesh) {
            UTOPIA_TRACE_REGION_BEGIN("FunctionSpace::init(mesh)");

            impl_->mesh = mesh;
            switch (mesh->spatial_dimension()) {
                case 2: {
                    using DMesh = ::mars::DistributedMesh<::mars::ElementType::Quad4>;
                    auto mesh_impl = mesh->raw_type<DMesh>();
                    impl_->init(*mesh_impl);
                    break;
                }
                case 3: {
                    using DMesh = ::mars::DistributedMesh<::mars::ElementType::Hex8>;
                    auto mesh_impl = mesh->raw_type<DMesh>();
                    impl_->init(*mesh_impl);
                    break;
                }
                default: {
                    assert(false);
                    Utopia::Abort("Trying to create mesh with unsupported dimension!");
                }
            }

            UTOPIA_TRACE_REGION_END("FunctionSpace::init(mesh)");
        }

        void FunctionSpace::update(const SimulationTime<Scalar> &time) { impl_->dirichlet_boundary.update(time); }

        bool FunctionSpace::write(const Path &path, const Vector &x) { return impl_->write(path, x); }

        void FunctionSpace::read(Input &in) {
            if (!impl_->mesh) {
                impl_->mesh = std::make_shared<Mesh>();
                in.get("mesh", *impl_->mesh);

                if (impl_->mesh->empty()) {
                    // Unit square
                    impl_->mesh->unit_cube(2, 2, 0);
                }
            }

            impl_->read_meta(in);
            init(impl_->mesh);
        }
        void FunctionSpace::describe(std::ostream &os) const {
            if (impl_->handler) {
                handler()->describe();
            } else {
                os << "FunctionSpace::describe: null\n";
            }
        }

        std::shared_ptr<FunctionSpace::Mesh> FunctionSpace::mesh_ptr() const { return impl_->mesh; }
        const FunctionSpace::Mesh &FunctionSpace::mesh() const { return *impl_->mesh; }
        FunctionSpace::Mesh &FunctionSpace::mesh() { return *impl_->mesh; }

        FunctionSpace::SizeType FunctionSpace::n_dofs() const { return handler()->n_dofs(); }
        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const { return handler()->n_local_dofs(); }
        int FunctionSpace::n_var() const { return impl_->n_var; }
        void FunctionSpace::set_n_var(const int n_var) { impl_->n_var = n_var; }

        void FunctionSpace::create_vector(Vector &v) const { v.zeros(layout(comm(), n_local_dofs(), n_dofs())); }
        void FunctionSpace::create_local_vector(Vector &) const { Utopia::Abort("IMPLEMENT ME"); }

        void FunctionSpace::create_matrix(Matrix &m) const {
            UTOPIA_TRACE_REGION_BEGIN("FunctionSpace::create_matrix(m)");

            using MapType = Matrix::MapType;

            handler()->ensure_sparsity_pattern();

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
                Teuchos::rcp(new Matrix::CrsMatrixType(handler()->new_crs_matrix(), row_map, col_map));
            m.wrap(mat, true);

            UTOPIA_TRACE_REGION_END("FunctionSpace::create_matrix(m)");
        }

        void FunctionSpace::copy_at_constrained_nodes(const Vector &in, Vector &out) const {
            UTOPIA_TRACE_REGION_BEGIN("FunctionSpace::copy_at_constrained_nodes(in, out)");

            for (auto &bc : impl_->dirichlet_boundary) {
                handler()->copy_at_constrained_nodes(in, out, bc->name, bc->component);
            }

            UTOPIA_TRACE_REGION_END("FunctionSpace::copy_at_constrained_nodes(in, out)");
        }

        void FunctionSpace::apply_constraints(Matrix &m, const Scalar diag_value) const {
            UTOPIA_TRACE_REGION_BEGIN("FunctionSpace::apply_constraints(m, diag_value)");
            for (auto &bc : impl_->dirichlet_boundary) {
                handler()->matrix_apply_constraints(m, diag_value, bc->name, bc->component);
            }
            UTOPIA_TRACE_REGION_END("FunctionSpace::apply_constraints(m, diag_value)");
        }

        void FunctionSpace::apply_constraints(Vector &v) const {
            UTOPIA_TRACE_REGION_BEGIN("FunctionSpace::apply_constraints(v)");

            for (auto &bc : impl_->dirichlet_boundary) {
                auto u_bc = std::dynamic_pointer_cast<DirichletBoundary::UniformCondition>(bc);

                if (u_bc) {
                    handler()->vector_apply_constraints(v, u_bc->value(), u_bc->name, u_bc->component);
                } else {
                    assert(false);
                }
            }

            UTOPIA_TRACE_REGION_END("FunctionSpace::apply_constraints(v)");
        }

        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) const {
            UTOPIA_TRACE_REGION_BEGIN("FunctionSpace::apply_constraints(m, v)");

            for (auto &bc : impl_->dirichlet_boundary) {
                auto u_bc = std::dynamic_pointer_cast<DirichletBoundary::UniformCondition>(bc);

                if (u_bc) {
                    handler()->matrix_apply_constraints(m, 1.0, u_bc->name, u_bc->component);
                    handler()->vector_apply_constraints(v, u_bc->value(), u_bc->name, u_bc->component);
                } else {
                    assert(false);
                }
            }

            UTOPIA_TRACE_REGION_END("FunctionSpace::apply_constraints(m, v)");
        }

        void FunctionSpace::apply_zero_constraints(Vector &vec) const {
            UTOPIA_TRACE_REGION_BEGIN("FunctionSpace::apply_zero_constraints(vec)");
            for (auto &bc : impl_->dirichlet_boundary) {
                handler()->vector_apply_constraints(vec, 0., bc->name, bc->component);
            }
            UTOPIA_TRACE_REGION_END("FunctionSpace::apply_zero_constraints(vec)");
        }

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name,
                                                             const Scalar &value,
                                                             const int component) {
            assert(component < n_var());
            DirichletBoundary::UniformCondition dirichlet_boundary{name, value, component};
            impl_->dirichlet_boundary.add(dirichlet_boundary);
        }

        bool FunctionSpace::empty() const { return !static_cast<bool>(impl_->handler); }

        void FunctionSpace::global_to_local(const Vector &global, Vector &local) const {}
        void FunctionSpace::local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const {}

        const std::string &FunctionSpace::name() const { return impl_->name; }

        const Factory &FunctionSpace::factory() const { return handler()->factory(); }

        template <class RawType>
        std::shared_ptr<RawType> FunctionSpace::raw_type() const {
            return impl_->raw_type<RawType>();
        }

    }  // namespace mars

}  // namespace utopia

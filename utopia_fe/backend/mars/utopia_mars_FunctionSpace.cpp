#include "utopia_mars_FunctionSpace.hpp"
#include <memory>

#include "utopia_DirichletBoundary.hpp"
#include "utopia_FEVar.hpp"
#include "utopia_mars_FEHandler.hpp"

#include "mars_pvtu_writer.hpp"

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
                handler_impl->init(mesh_impl);
                handler = handler_impl;

#ifdef UTOPIA_WITH_VTK
                write = [handler_impl](const Path &path, const Vector &x) -> bool {
                    ::mars::PVTUMeshWriter<typename FEHandler::DofHandler, typename FEHandler::FEDofMap> w;
                    auto x_kokkos = x.raw_type()->getLocalViewHost();
                    return w.write_vtu(
                        path.to_string(), handler_impl->get_dof_handler(), handler_impl->get_fe_dof_map(), x_kokkos);
                };
#else
                write = (const Path &, const Vector &)->bool { return false; };
#endif  // UTOPIA_WITH_VTK
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
            int n_var{0};
        };

        std::shared_ptr<IFEHandler> FunctionSpace::handler() const { return impl_->handler; }

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
        }

        bool FunctionSpace::write(const Path &path, const Vector &x) { return impl_->write(path, x); }

        void FunctionSpace::read(Input &in) {
            if (impl_->mesh) {
                impl_->mesh = std::make_shared<Mesh>();
                in.get("mesh", *impl_->mesh);
            }

            impl_->read_meta(in);
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
                Teuchos::rcp(new Matrix::CrsMatrixType(handler()->new_crs_matrix(), row_map, col_map));
            m.wrap(mat, true);
        }

        void FunctionSpace::apply_constraints(Matrix &m, const Scalar diag_value) {
            for (auto &bc : impl_->dirichlet_boundary.conditions) {
                handler()->matrix_apply_constraints(m, diag_value, bc.name);
            }
        }

        void FunctionSpace::apply_constraints(Vector &v) {
            for (auto &bc : impl_->dirichlet_boundary.conditions) {
                handler()->vector_apply_constraints(v, bc.value, bc.name);
            }
        }

        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) {
            for (auto &bc : impl_->dirichlet_boundary.conditions) {
                handler()->matrix_apply_constraints(m, 1.0, bc.name);
                handler()->vector_apply_constraints(v, bc.value, bc.name);
            }
        }

        void FunctionSpace::apply_zero_constraints(Vector &vec) const {
            for (auto &bc : impl_->dirichlet_boundary.conditions) {
                handler()->vector_apply_constraints(vec, bc.value, bc.name);
            }
        }

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name,
                                                             const Scalar &value,
                                                             const int component) {
            assert(component < n_var());
            DirichletBoundary::Condition dirichlet_boundary{name, value, component};
            impl_->dirichlet_boundary.conditions.push_back(dirichlet_boundary);
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
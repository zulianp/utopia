#include "utopia_mars_FunctionSpace.hpp"

#include "utopia_DirichletBoundary.hpp"
#include "utopia_FEVar.hpp"

#include "mars.hpp"
#include "mars_base.hpp"
#include "mars_context.hpp"
#include "mars_execution_context.hpp"

namespace utopia {
    namespace mars {

        class FunctionSpace::Impl {
        public:
            using DeviceType = typename Kokkos::Device<Kokkos::DefaultExecutionSpace, KokkosSpace>;
            // using Ordinal = unsigned long;
            using Ordinal = SizeType;

            using MarsCrsMatrix = typename KokkosSparse::CrsMatrix<Scalar, Ordinal, DeviceType, void, Ordinal>;

            template <class DMesh>
            void init(DMesh &mesh_impl) {
                static constexpr ::mars::Integer Degree = 1;
                using DofHandler = ::mars::DofHandler<DMesh, Degree>;
                using FEDofMap = ::mars::FEDofMap<DofHandler>;
                // using DM = ::mars::DM<DofHandler, Scalar, Scalar, Scalar>;
                using FE = ::mars::FEDofMap<DofHandler>;
                using SPattern = ::mars::SparsityPattern<Scalar, ::mars::Integer, Ordinal, DofHandler>;

                auto dof_handler_impl = std::make_shared<DofHandler>(&mesh_impl, mesh->raw_type_context());
                dof_handler_impl->enumerate_dofs();

                dof_handler = dof_handler_impl;

                auto fe_dof_map_impl = std::make_shared<FEDofMap>(build_fe_dof_map(*dof_handler_impl));
                fe_dof_map = fe_dof_map_impl;

                SPattern sp(*dof_handler_impl);
                sp.build_pattern(*fe_dof_map_impl);
                crs_matrix = sp.get_matrix();

                describe = [dof_handler_impl, fe_dof_map_impl]() {
                    dof_handler_impl->print_dofs();
                    // ::mars::print_fe_dof_map(*dof_handler_impl, *fe_dof_map_impl);
                };
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

            std::shared_ptr<::mars::IDofHandler> dof_handler;
            std::shared_ptr<::mars::IFEDofMap> fe_dof_map;
            std::shared_ptr<Mesh> mesh;

            std::string name;
            DirichletBoundary dirichlet_boundary;
            std::vector<FEVar> variables;

            MarsCrsMatrix crs_matrix;

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

        std::shared_ptr<FunctionSpace::Mesh> FunctionSpace::mesh_ptr() const {}
        const FunctionSpace::Mesh &FunctionSpace::mesh() const {}
        FunctionSpace::Mesh &FunctionSpace::mesh() {}

        FunctionSpace::SizeType FunctionSpace::n_dofs() const {}
        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const {}
        int FunctionSpace::n_var() const {}
        void FunctionSpace::set_n_var(const int n_var) {}

        void FunctionSpace::create_vector(Vector &v) const {}
        void FunctionSpace::create_local_vector(Vector &v) const {}

        void FunctionSpace::create_matrix(Matrix &m) const {
            Matrix::RCPCrsMatrixType mat;
            // mat.reset(new Matrix::CrsMatrixType());

            m.wrap(mat, true);
        }

        void FunctionSpace::apply_constraints(Matrix &m, const Scalar diag_value) {}
        void FunctionSpace::apply_constraints(Vector &v) {}
        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) {}
        void FunctionSpace::apply_zero_constraints(Vector &vec) const {}

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name,
                                                             const Scalar &value,
                                                             const int component) {}

        bool FunctionSpace::empty() const {}

        void FunctionSpace::global_to_local(const Vector &global, Vector &local) const {}
        void FunctionSpace::local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const {}

        const std::string &FunctionSpace::name() const {}

    }  // namespace mars

}  // namespace utopia

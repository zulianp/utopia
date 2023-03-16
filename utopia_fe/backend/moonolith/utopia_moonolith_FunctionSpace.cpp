#include "utopia_moonolith_FunctionSpace.hpp"

#include "moonolith_function_space.hpp"
#include "moonolith_mesh.hpp"
#include "moonolith_vtu_mesh_writer.hpp"

#include <memory>

namespace utopia {
    namespace moonolith {
        template <int Dim>
        using MoonolithMesh = ::moonolith::Mesh<FunctionSpace::Scalar, Dim>;

        template <int Dim>
        using MoonolithFunctionSpace = ::moonolith::FunctionSpace<MoonolithMesh<Dim>>;

        class FunctionSpace::Impl {
        public:
            template <int Dim>
            void wrap(const std::shared_ptr<MoonolithFunctionSpace<Dim>> &ptr) {
                space = ptr;

                manifold_dim = [ptr]() -> int { return Dim; };
                n_dofs = [ptr]() -> SizeType { return ptr->dof_map().n_dofs(); };
                n_local_dofs = [ptr]() -> SizeType { return ptr->dof_map().n_local_dofs(); };
                write = [ptr](const Path &path, const Vector &x) -> bool {
                    ::moonolith::VTUMeshWriter<MoonolithMesh<Dim>> writer;
                    auto x_view = local_view_device(x);
                    return writer.write(path.to_string(), ptr->mesh(), x_view);
                };

                describe = [ptr](std::ostream &os) { os << "dim: " << Dim; };
            }

            std::shared_ptr<Mesh> mesh;
            std::shared_ptr<::moonolith::IFunctionSpace> space;

            // Prototype based programming
            std::function<int()> manifold_dim;
            std::function<SizeType()> n_dofs;
            std::function<SizeType()> n_local_dofs;
            std::function<bool(const Path &, const Vector &)> write;
            std::function<void(std::ostream &os)> describe;
        };

        FunctionSpace::FunctionSpace(const Comm &comm) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = std::make_shared<Mesh>(comm);
        }

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = mesh;
        }

        FunctionSpace::~FunctionSpace() = default;

        bool FunctionSpace::write(const Path &path, const Vector &x) { return impl_->write(path, x); }

        void FunctionSpace::init(const std::shared_ptr<Mesh> &mesh, const bool init_as_iso_paramatric) {
            impl_->mesh = mesh;
            switch (mesh->spatial_dimension()) {
                case 1: {
                    auto space = std::make_shared<MoonolithFunctionSpace<1>>(mesh->raw_type<1>());

                    if (init_as_iso_paramatric) {
                        space->make_iso_parametric();
                    }

                    space->dof_map().set_n_local_dofs(mesh->n_local_nodes());
                    space->dof_map().set_n_dofs(mesh->n_nodes());
                    impl_->wrap(space);
                    break;
                }

                case 2: {
                    auto space = std::make_shared<MoonolithFunctionSpace<2>>(mesh->raw_type<2>());

                    if (init_as_iso_paramatric) {
                        space->make_iso_parametric();
                    }

                    space->dof_map().set_n_local_dofs(mesh->n_local_nodes());
                    space->dof_map().set_n_dofs(mesh->n_nodes());
                    impl_->wrap(space);
                    break;
                }

                case 3: {
                    auto space = std::make_shared<MoonolithFunctionSpace<3>>(mesh->raw_type<3>());

                    if (init_as_iso_paramatric) {
                        space->make_iso_parametric();
                    }

                    space->dof_map().set_n_local_dofs(mesh->n_local_nodes());
                    space->dof_map().set_n_dofs(mesh->n_nodes());
                    impl_->wrap(space);
                    break;
                }

                default: {
                    utopia::err() << "FunctionSpace::read: unsupported dimension\n";
                    Utopia::Abort();
                    break;
                }
            }
        }

        void FunctionSpace::read(Input &in) {
            auto mesh = std::make_shared<Mesh>();
            mesh->read(in);
            init(mesh);
        }

        void FunctionSpace::create_vector(Vector &v) const { v.zeros(layout(comm(), n_local_dofs(), n_dofs())); }

        void FunctionSpace::describe(std::ostream &os) const { impl_->describe(os); }

        std::shared_ptr<Mesh> FunctionSpace::mesh_ptr() const { return impl_->mesh; }

        const Mesh &FunctionSpace::mesh() const {
            assert(impl_->mesh);
            return *impl_->mesh;
        }

        Mesh &FunctionSpace::mesh() {
            assert(impl_->mesh);
            return *impl_->mesh;
        }

        template <int Dim>
        std::shared_ptr<::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, Dim>>>
        FunctionSpace::raw_type() const {
            // assert(Dim == impl_->manifold_dim());
            return std::dynamic_pointer_cast<MoonolithFunctionSpace<Dim>>(impl_->space);
        }

        template <int Dim>
        void FunctionSpace::wrap(
            const std::shared_ptr<::moonolith::FunctionSpace<::moonolith::Mesh<Scalar, Dim>>> &space) {
            impl_->mesh->wrap(space->mesh_ptr());
            impl_->wrap(space);
        }

        FunctionSpace::SizeType FunctionSpace::n_dofs() const { return impl_->n_dofs(); }
        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const { return impl_->n_local_dofs(); }

        // Explicit instantiations
        template std::shared_ptr<::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, 1>>>
        FunctionSpace::raw_type() const;

        template std::shared_ptr<::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, 2>>>
        FunctionSpace::raw_type() const;

        template std::shared_ptr<::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, 3>>>
        FunctionSpace::raw_type() const;

        template void FunctionSpace::wrap(
            const std::shared_ptr<::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, 1>>> &space);

        template void FunctionSpace::wrap(
            const std::shared_ptr<::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, 2>>> &space);

        template void FunctionSpace::wrap(
            const std::shared_ptr<::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, 3>>> &space);

    }  // namespace moonolith
}  // namespace utopia

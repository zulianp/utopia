#ifndef UTOPIA_SDF_OBSTACLE_IMPL_HPP
#define UTOPIA_SDF_OBSTACLE_IMPL_HPP

#include "utopia.hpp"
#include "utopia_IPTransfer.hpp"
#include "utopia_SymbolicFunction.hpp"

#include "utopia_Field.hpp"
#include "utopia_SDFObstacle.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_moonolith_HouseholderReflection.hpp"

#include "utopia_sfem_Mesh.hpp"
#include "utopia_sfem_SDF.hpp"

namespace utopia {

    template <class FunctionSpace>
    class SDFObstacle<FunctionSpace>::Impl {
    public:
        std::shared_ptr<Field> gap;
        std::shared_ptr<GradientField> normals;
        Matrix orthogonal_trafo;
        Vector is_contact;
        Scalar infinity{1e10};
        utopia::sfem::SDF sdf;
        utopia::sfem::Mesh mesh;

        void resample_to_mesh_surface_of(FunctionSpace &space) {
            typename FunctionSpace::Mesh surface_mesh;
            extract_surface(space.mesh(), mesh);
            Vector surface_gap, surface_normals;
            sdf.to_mesh(mesh, surface_gap, surface_normals);

            // Convert to space context
            auto gap = std::make_shared<Field>();
            auto normals = std::make_shared<GradientField>();

            space.create_field(*gap);
            space.create_field(*normals);

            gap->data().set(infinity);
            space.create_vector(is_contact);

            {
                // Scope for views
                auto surface_gap_view = const_local_view_device(surface_gap);
                auto gap_view = local_view_device(gap->data());

                auto surface_normals_view = const_local_view_device(surface_normals);
                auto normals_view = local_view_device(normals->data());

                auto is_contact_view = local_view_device(is_contact);

                auto node_mapping = mesh.node_mapping();

                const int n_var = space.n_var();
                for (ptrdiff_t i = 0; i < mesh.n_local_nodes(); i++) {
                    SizeType node = node_mapping[i];

                    gap_view.set(node * n_var, surface_gap_view.get(i));
                    normals_view.set(node * n_var, surface_normals_view.get(i * 3));
                    normals_view.set(node * n_var + 1, surface_normals_view.get(i * 3 + 1));
                    normals_view.set(node * n_var + 2, surface_normals_view.get(i * 3 + 2));

                    is_contact_view.set(node * n_var, 1);
                }
            }

            space.apply_zero_constraints(is_contact);
            HouseholderReflectionForContact<Matrix, Vector, 3>::build(is_contact, normals->data(), orthogonal_trafo);
        }
    };

    template <class FunctionSpace>
    SDFObstacle<FunctionSpace>::SDFObstacle() : impl_(utopia::make_unique<Impl>()) {}

    template <class FunctionSpace>
    SDFObstacle<FunctionSpace>::~SDFObstacle() = default;

    template <class FunctionSpace>
    void SDFObstacle<FunctionSpace>::read(Input &in) {
        impl_->sdf.read(in);
        // in.get("field_rescale", impl_->field_rescale);
        // in.get("field_offset", impl_->field_offset);
    }

    template <class FunctionSpace>
    void SDFObstacle<FunctionSpace>::describe(std::ostream &os) const {
        UTOPIA_UNUSED(os);
    }

    template <class FunctionSpace>
    bool SDFObstacle<FunctionSpace>::assemble(FunctionSpace &space) {
        impl_->resample_to_mesh_surface_of(space);
        return true;
    }

    template <class FunctionSpace>
    void SDFObstacle<FunctionSpace>::transform(const Matrix &in, Matrix &out) {
        out = transpose(impl_->orthogonal_trafo) * in * impl_->orthogonal_trafo;
    }

    template <class FunctionSpace>
    void SDFObstacle<FunctionSpace>::transform(const Vector &in, Vector &out) {
        out = transpose(impl_->orthogonal_trafo) * in;
    }

    template <class FunctionSpace>
    void SDFObstacle<FunctionSpace>::inverse_transform(const Vector &in, Vector &out) {
        out = impl_->orthogonal_trafo * in;
    }

    template <class FunctionSpace>
    const typename SDFObstacle<FunctionSpace>::Vector &SDFObstacle<FunctionSpace>::gap() const {
        assert(impl_->gap);
        return impl_->gap->data();
    }

    template <class FunctionSpace>
    const typename SDFObstacle<FunctionSpace>::Vector &SDFObstacle<FunctionSpace>::is_contact() const {
        return impl_->is_contact;
    }

    template <class FunctionSpace>
    const typename SDFObstacle<FunctionSpace>::Vector &SDFObstacle<FunctionSpace>::normals() const {
        assert(impl_->normals);
        return impl_->normals->data();
    }

    template <class FunctionSpace>
    std::shared_ptr<typename Traits<FunctionSpace>::Matrix> SDFObstacle<FunctionSpace>::orthogonal_transformation() {
        return make_ref(impl_->orthogonal_trafo);
    }

}  // namespace utopia

#endif  // UTOPIA_SDF_OBSTACLE_IMPL_HPP

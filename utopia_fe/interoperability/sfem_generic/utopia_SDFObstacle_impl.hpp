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

        void resample_to_mesh_surface_of(const FunctionSpace &space) {
            typename FunctionSpace::Mesh surface_mesh;
            extract_surface(space.mesh(), mesh);
            Vector surface_gap, surface_normal;
            sdf.to_mesh(mesh, surface_gap, surface_normal);

            // surface to volume dof mapping
            // TODO
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
        // auto gap = std::make_shared<Field>();
        // auto normals = std::make_shared<GradientField>();
        // int dim = space.mesh().spatial_dimension();

        // space.create_field(*gap);
        // space.create_field(*normals);

        // assert(!gap->empty());
        // assert(!normals->empty());

        // impl_->gap = gap;
        // impl_->normals = normals;
        // impl_->is_contact.zeros(layout(normals->data()));

        // impl_->gap->data().set(impl_->infinity);

        // int n_var = space.n_var();

        // {
        //     auto g_view = view_device(impl_->gap->data());
        //     auto n_view = view_device(impl_->normals->data());
        //     auto c_view = view_device(impl_->is_contact);

        //     Range r = range(impl_->gap->data());
        //     SizeType r_begin = r.begin() / n_var;
        //     SizeType r_end = r.end() / n_var;

        //     for (auto &f : impl_->functions) {
        //         auto fun = [this, dim, g_view, n_view, c_view, n_var, r_begin, r_end, &f](const SizeType idx,
        //                                                                                   const Scalar *point) {
        //             // if (idx < r_begin || idx >= r_end) return;

        //             Scalar point3[3] = {0.0, 0.0, 0.0};
        //             Scalar normal3[3] = {0.0, 0.0, 0.0};
        //             Scalar g = 0;

        //             for (int d = 0; d < dim; ++d) {
        //                 point3[d] = point[d];
        //             }

        //             f->eval(point3, g, normal3);
        //             g_view.set(idx * n_var, g);
        //             c_view.set(idx * n_var, 1);

        //             for (int d = 0; d < dim; ++d) {
        //                 n_view.set(idx * n_var + d, normal3[d]);
        //             }
        //         };

        //         if (f->side.empty()) {
        //             space.node_eval(fun);
        //         } else {
        //             space.node_eval(f->side, fun);
        //         }
        //     }
        // }

        // space.apply_zero_constraints(impl_->is_contact);
        // space.apply_zero_constraints(normals->data());

        // switch (dim) {
        //     // case 2: {
        //     //     HouseholderReflectionForContact<Matrix, Vector, 2>::build(
        //     //         impl_->is_contact, normals->data(), impl_->orthogonal_trafo);
        //     //     break;
        //     // }
        //     case 3: {
        //         HouseholderReflectionForContact<Matrix, Vector, 3>::build(
        //             impl_->is_contact, normals->data(), impl_->orthogonal_trafo);
        //         break;
        //     }
        //     default: {
        //         Utopia::Abort("Unsupported dimension!");
        //         return false;
        //     }
        // }

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

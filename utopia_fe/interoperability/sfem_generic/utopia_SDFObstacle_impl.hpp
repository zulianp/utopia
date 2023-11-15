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

#include "utopia_fe_Core.hpp"

namespace utopia {

    template <class FunctionSpace>
    class SDFObstacle<FunctionSpace>::Impl {
    public:
        std::shared_ptr<Field> gap;
        std::shared_ptr<GradientField> normals;
        Matrix orthogonal_trafo;
        Vector is_contact;
        Scalar infinity{1e10};
        Scalar shift{0};
        Scalar cutoff{1e8};
        utopia::sfem::SDF sdf;
        utopia::sfem::Mesh mesh;
        bool export_gap{false};
        bool verbose{false};

        void normalize(Vector &normal) {
            ReadAndWrite<Vector> rw_normal(normal);
            auto r = normal.range();

            const int dim = mesh.spatial_dimension();
            Scalar n[3];
            for (auto i = r.begin(); i < r.end(); i += dim) {
                Scalar len_n = 0;
                for (int d = 0; d < dim; ++d) {
                    n[d] = normal.get(i + d);
                    len_n += n[d] * n[d];
                }

                len_n = std::sqrt(len_n);

                if (len_n == 0.0) continue;

                for (int d = 0; d < dim; ++d) {
                    n[d] /= len_n;
                }

                for (int d = 0; d < dim; ++d) {
                    normal.set(i + d, n[d]);
                }
            }
        }

        void resample_to_mesh_surface_of(FunctionSpace &space) {
            typename FunctionSpace::Mesh surface_mesh;
            extract_surface(space.mesh(), mesh);

            if (0) {
                static int export_mesh = 0;
                mesh.write("dbg_mesh" + std::to_string(export_mesh++));
            }

            Vector surface_gap, surface_normals;
            sdf.to_mesh(mesh, surface_gap, surface_normals);
            surface_gap.shift(shift);

            Scalar norm_gap = norm2(surface_gap);

            if (!space.comm().rank()) {
                utopia::out() << "Surface Mesh #nodes " << mesh.n_local_nodes() << "\n";
                utopia::out() << "Norm Gap: " << norm_gap << "\n";
            }

            // Convert to space context
            gap = std::make_shared<Field>();
            normals = std::make_shared<GradientField>();

            space.create_field(*gap);
            space.create_field(*normals);

            gap->data().set(infinity);
            space.create_vector(is_contact);

            auto &&local_to_global = space.dof_map().local_to_global();

            if (sdf.has_weights()) {
                Vector weights;

                Vector local_gap;
                Vector local_normals;
                Vector local_is_contact;
                Vector local_weights;

                space.create_local_vector(local_gap);
                space.create_local_vector(local_normals);
                space.create_local_vector(local_is_contact);

                local_gap.set(infinity);
                local_is_contact.set(0);
                local_normals.set(0);

                // Scope for views
                {
                    auto surface_gap_view = const_local_view_device(surface_gap);
                    auto surface_normals_view = const_local_view_device(surface_normals);
                    auto surface_weights_view = local_view_device(sdf.weights());

                    auto gap_view = local_view_device(local_gap);
                    auto normals_view = local_view_device(local_normals);
                    auto weights_view = local_view_device(local_weights);
                    auto is_contact_view = local_view_device(local_is_contact);

                    auto node_mapping = mesh.node_mapping();

                    auto rr = gap->data().range();
                    const int n_var = space.n_var();
                    const SizeType local_size = gap->data().local_size();
                    for (ptrdiff_t i = 0; i < mesh.n_local_nodes(); i++) {
                        const Scalar gg = surface_gap_view.get(i);
                        if (gg > cutoff) continue;

                        SizeType node = node_mapping[i];
                        // Offset for dof number
                        node *= n_var;

                        assert(node + 2 < is_contact.local_size());
                        assert(node + 2 < gap->data().local_size());
                        assert(node + 2 < normals->data().local_size());

                        gap_view.set(node, gg);
                        weights_view.set(node, surface_weights_view.get(i));
                        normals_view.set(node, surface_normals_view.get(i * 3));
                        normals_view.set(node + 1, surface_normals_view.get(i * 3 + 1));
                        normals_view.set(node + 2, surface_normals_view.get(i * 3 + 2));
                        is_contact_view.set(node, 1);
                    }
                }

                space.local_to_global(local_gap, gap->data(), utopia::ADD_MODE);
                space.local_to_global(local_normals, normals->data(), utopia::ADD_MODE);
                space.local_to_global(local_is_contact, is_contact, utopia::ADD_MODE);
                space.local_to_global(local_weights, weights, utopia::ADD_MODE);

                // Clamp
                is_contact.e_min(1);

                // Divide by weight
                e_pseudo_inv(gap->data(), weights, 1e-16);
                normalize(normals->data());

            } else if (!local_to_global.empty()) {
                // Scope for views
                {
                    auto surface_gap_view = const_local_view_device(surface_gap);
                    auto gap_view = local_view_device(gap->data());

                    auto surface_normals_view = const_local_view_device(surface_normals);
                    auto normals_view = local_view_device(normals->data());

                    auto is_contact_view = local_view_device(is_contact);

                    auto node_mapping = mesh.node_mapping();

                    auto rr = gap->data().range();
                    const int n_var = space.n_var();
                    const SizeType local_size = gap->data().local_size();
                    for (ptrdiff_t i = 0; i < mesh.n_local_nodes(); i++) {
                        const Scalar gg = surface_gap_view.get(i);
                        if (gg > cutoff) continue;

                        SizeType node = node_mapping[i];
                        node = local_to_global(node, 0) - rr.begin();

                        if (node >= local_size || node < 0) {
                            // Skip ghost nodes!
                            continue;
                        }

                        assert(node + 2 < is_contact.local_size());
                        assert(node + 2 < gap->data().local_size());
                        assert(node + 2 < normals->data().local_size());

                        gap_view.set(node, gg);
                        normals_view.set(node, surface_normals_view.get(i * 3));
                        normals_view.set(node + 1, surface_normals_view.get(i * 3 + 1));
                        normals_view.set(node + 2, surface_normals_view.get(i * 3 + 2));
                        is_contact_view.set(node, 1);
                    }
                }

            } else {
                // Serial case
                auto surface_gap_view = const_local_view_device(surface_gap);
                auto gap_view = local_view_device(gap->data());

                auto surface_normals_view = const_view_device(surface_normals);
                auto normals_view = view_device(normals->data());

                auto is_contact_view = view_device(is_contact);

                auto node_mapping = mesh.node_mapping();

                const int n_var = space.n_var();
                for (ptrdiff_t i = 0; i < mesh.n_local_nodes(); i++) {
                    const Scalar gg = surface_gap_view.get(i);
                    if (gg > cutoff) continue;

                    SizeType node = node_mapping[i];
                    // Offset for dof number
                    node *= n_var;

                    gap_view.set(node, gg);
                    normals_view.set(node, surface_normals_view.get(i * 3));
                    normals_view.set(node + 1, surface_normals_view.get(i * 3 + 1));
                    normals_view.set(node + 2, surface_normals_view.get(i * 3 + 2));
                    is_contact_view.set(node, 1);
                }
            }

            if (verbose) {
                SizeType num_contacts = sum(is_contact);
                if (!space.comm().rank()) {
                    utopia::out() << "SDFObstacle: num_contacts = " << num_contacts << "\n";
                }
            }

            space.apply_zero_constraints(is_contact);
            HouseholderReflectionForContact<Matrix, Vector, 3>::build(is_contact, normals->data(), orthogonal_trafo);

            if (export_gap) {
                static int count = 0;

                Vector director = normals->data();

                {
                    // Scope for views
                    auto gap_view = local_view_device(gap->data());
                    auto director_view = local_view_device(director);
                    auto is_contact_view = local_view_device(is_contact);

                    parallel_for(
                        local_range_device(director), UTOPIA_LAMBDA(const SizeType i) {
                            const SizeType ii = (i / 3) * 3;
                            const Scalar g = gap_view.get(ii);
                            // const Scalar g = 1;
                            const Scalar ind = is_contact_view.get(ii);

                            director_view.set(i, director_view.get(i) * g * ind);
                        });
                }

                space.write("director_" + std::to_string(count++) + ".e", director);
            }
        }
    };

    template <class FunctionSpace>
    SDFObstacle<FunctionSpace>::SDFObstacle() : impl_(utopia::make_unique<Impl>()) {}

    template <class FunctionSpace>
    SDFObstacle<FunctionSpace>::~SDFObstacle() = default;

    template <class FunctionSpace>
    void SDFObstacle<FunctionSpace>::read(Input &in) {
        impl_->sdf.read(in);
        in.get("infinity", impl_->infinity);
        in.get("export_gap", impl_->export_gap);
        in.get("shift", impl_->shift);
        in.get("cutoff", impl_->cutoff);
        in.get("verbose", impl_->verbose);
    }

    template <class FunctionSpace>
    void SDFObstacle<FunctionSpace>::describe(std::ostream &os) const {
        UTOPIA_UNUSED(os);
    }

    template <class FunctionSpace>
    bool SDFObstacle<FunctionSpace>::assemble(FunctionSpace &space) {
        UTOPIA_TRACE_SCOPE("SDFObstacle::assemble");
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

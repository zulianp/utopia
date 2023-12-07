#ifndef UTOPIA_SDF_OBSTACLE_IMPL_HPP
#define UTOPIA_SDF_OBSTACLE_IMPL_HPP

#include <utopia_DeviceView.hpp>
#include <utopia_Enums.hpp>
#include <utopia_Epsilon.hpp>
#include <utopia_RangeDevice.hpp>
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
        std::vector<std::string> exclude;

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
            extract_surface(space.mesh(), mesh, exclude);

            if (0) {
                static int export_mesh = 0;
                mesh.write("dbg_mesh" + std::to_string(export_mesh++));
            }

            Vector surface_gap, surface_normals;
            sdf.to_mesh(mesh, surface_gap, surface_normals);

            if (verbose) {
                Scalar sum_gap = sum(surface_gap);
                Scalar max_gap = max(surface_gap);
                Scalar min_gap = min(surface_gap);

                std::stringstream ss;
                ss << "Surface Mesh #nodes " << mesh.n_local_nodes() << "\n";
                ss << "sum_gap: " << sum_gap << "\n";
                ss << "min_gap: " << min_gap << "\n";
                ss << "max_gap: " << max_gap << "\n";
                ss << "shift: " << shift << "\n";
                space.comm().root_print(ss.str());
            }

            // Convert to space context
            gap = std::make_shared<Field>();
            normals = std::make_shared<GradientField>();

            space.create_field(*gap);
            space.create_field(*normals);
            space.create_vector(is_contact);
            is_contact.set(0);
            normals->data().set(0);

            auto &&local_to_global = space.dof_map().local_to_global();

            if (sdf.has_weights()) {
                UTOPIA_TRACE_SCOPE("SDFObstacle::resample_to_mesh_surface_of: parallel/projection");

                Vector weights;
                space.create_vector(weights);
                weights.set(0.);
                gap->data().set(0.);

                // Scope for views
                {
                    auto surface_gap_view = const_local_view_host(surface_gap);
                    auto surface_normals_view = const_local_view_host(surface_normals);
                    auto surface_weights_view = const_local_view_host(sdf.weights());
                    auto node_mapping = mesh.node_mapping();

                    Write<Vector>                                   //
                        w_g(gap->data(), utopia::GLOBAL_ADD),       //
                        w_sn(normals->data(), utopia::GLOBAL_ADD),  //
                        w_w(weights, utopia::GLOBAL_ADD),           //
                        w_ic(is_contact, utopia::GLOBAL_ADD);

                    auto &g = gap->data();
                    auto &n = normals->data();

                    const bool has_l2g = !local_to_global.empty();
                    const int n_var = space.n_var();

                    auto rr = gap->data().range();
                    for (SizeType i = 0; i < surface_gap.size(); i++) {
                        const Scalar value = surface_gap_view.get(i);
                        const SizeType node = node_mapping[i];
                        SizeType globalid;

                        if (has_l2g) {
                            globalid = local_to_global(node, 0);
                        } else {
                            globalid = node * n_var;
                        }

                        g.c_add(globalid, value);
                        n.c_add(globalid + 0, surface_normals_view.get(i * 3 + 0));
                        n.c_add(globalid + 1, surface_normals_view.get(i * 3 + 1));
                        n.c_add(globalid + 2, surface_normals_view.get(i * 3 + 2));
                        weights.c_add(globalid, surface_weights_view.get(i));
                        is_contact.c_add(globalid, 1);
                    }
                }

                // Divide by weight
                e_pseudo_inv(weights, weights, device::Epsilon<Scalar>::value());
                gap->data() = e_mul(gap->data(), weights);
                normalize(normals->data());

                // if (0)  //
                {
                    RangeDevice<Vector> block_range(0, gap->data().local_size() / 3);

                    auto gap_view = local_view_device(gap->data());
                    auto normals_view = local_view_device(normals->data());
                    auto is_contact_view = local_view_device(is_contact);

                    parallel_for(
                        block_range, UTOPIA_LAMBDA(SizeType block_id) {
                            auto i = block_id * 3;
                            auto gg = gap_view.get(i);
                            if (gg > cutoff) {
                                is_contact_view.set(i, 0);
                            }

                            if (is_contact_view.get(i) == 0) {
                                normals_view.set(i, 0);
                                normals_view.set(i + 1, 0);
                                normals_view.set(i + 2, 0);
                                gap_view.set(i, infinity);
                            } else {
                                // Clamp to 1
                                is_contact_view.set(i, 1);
                            }

                            gap_view.set(i + 1, infinity);
                            gap_view.set(i + 2, infinity);
                        });
                }

            } else if (!local_to_global.empty()) {
                UTOPIA_TRACE_SCOPE("SDFObstacle::resample_to_mesh_surface_of: parallel");

                gap->data().set(infinity);

                // Scope for views
                {
                    auto surface_gap_view = const_local_view_host(surface_gap);
                    auto surface_normals_view = const_local_view_host(surface_normals);
                    auto node_mapping = mesh.node_mapping();

                    Write<Vector>                                      //
                        w_g(gap->data(), utopia::GLOBAL_INSERT),       //
                        w_sn(normals->data(), utopia::GLOBAL_INSERT),  //
                        w_ic(is_contact, utopia::GLOBAL_INSERT);

                    auto rr = gap->data().range();
                    for (ptrdiff_t i = 0; i < mesh.n_local_nodes(); i++) {
                        const Scalar gg = surface_gap_view.get(i);
                        if (gg > cutoff) continue;

                        SizeType globalid = local_to_global(node_mapping[i], 0);

                        gap->data().c_set(globalid, gg);
                        normals->data().c_set(globalid, surface_normals_view.get(i * 3));
                        normals->data().c_set(globalid + 1, surface_normals_view.get(i * 3 + 1));
                        normals->data().c_set(globalid + 2, surface_normals_view.get(i * 3 + 2));
                        is_contact.c_set(globalid, 1);
                    }
                }

            } else {
                UTOPIA_TRACE_SCOPE("SDFObstacle::resample_to_mesh_surface_of: serial");
                gap->data().set(infinity);

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

            gap->data().shift(-shift);

            space.apply_zero_constraints(is_contact);
            HouseholderReflectionForContact<Matrix, Vector, 3>::build(is_contact, normals->data(), orthogonal_trafo);

            if (export_gap) {
                static int count = 0;

                Vector director = normals->data();

                // if (0)  //
                {
                    // Scope for views
                    auto gap_view = local_view_device(gap->data());
                    auto director_view = local_view_device(director);
                    auto is_contact_view = local_view_device(is_contact);

                    parallel_for(
                        local_range_device(director), UTOPIA_LAMBDA(const SizeType i) {
                            const SizeType ii = (i / 3) * 3;
                            const Scalar g = gap_view.get(ii);
                            const Scalar ind = is_contact_view.get(ii);
                            director_view.set(i, director_view.get(i) * g * ind);
                        });
                }

                Vector filtered_gap = e_mul(is_contact, gap->data());
                // Vector filtered_gap = gap->data();

                space.write("gap_" + std::to_string(count) + ".e", filtered_gap);
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

        in.get("exclude", [&](Input &list) {
            list.get_all([&](Input &node) {
                std::string name;
                node.require("name", name);
                impl_->exclude.push_back(name);
            });
        });
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

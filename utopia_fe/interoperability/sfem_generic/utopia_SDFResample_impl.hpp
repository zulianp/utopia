#ifndef UTOPIA_SDF_RESAMPLE_IMPL_HPP
#define UTOPIA_SDF_RESAMPLE_IMPL_HPP

#include "utopia_Epsilon.hpp"
#include "utopia_Input.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_sfem_Mesh.hpp"
#include "utopia_sfem_SDF.hpp"

#include "utopia_SDFResample.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_stk_DofMap.hpp"

#include <vector>

namespace utopia {

    template <class FunctionSpace>
    class SDFResample<FunctionSpace>::Impl : public Configurable {
    public:
        utopia::sfem::SDF sdf;
        utopia::sfem::Mesh mesh;
        bool verbose{false};
        std::vector<std::string> exclude;

        Scalar shift{0};

        void read(Input &in) override {
            sdf.read(in);

            in.get("verbose", verbose);
            in.get("shift", shift);

            in.get("exclude", [&](Input &list) {
                list.get_all([&](Input &node) {
                    std::string name;
                    node.require("name", name);
                    exclude.push_back(name);
                });
            });
        }

        void resample_to_mesh_surface_of(FunctionSpace &space, Vector &out) {
            typename FunctionSpace::Mesh surface_mesh;
            extract_surface(space.mesh(), mesh, exclude);

            Vector surface_field;
            sdf.to_mesh(mesh, surface_field);
            surface_field.shift(shift);

            Scalar norm_sdf = norm2(surface_field);

            if (!space.comm().rank()) {
                utopia::out() << "Surface Mesh #nodes " << mesh.n_local_nodes() << "\n";
                utopia::out() << "Norm Gap: " << norm_sdf << "\n";
            }

            space.create_vector(out);

            auto &&local_to_global = space.dof_map().local_to_global();

            if (sdf.has_weights()) {
                Vector weights;

                Vector local_field;
                Vector local_weights;

                space.create_local_vector(local_field);
                space.create_local_vector(local_weights);

                local_field.set(0);
                local_weights.set(0);

                // Scope for views
                {
                    auto surface_field_view = const_local_view_host(surface_field);
                    auto surface_weights_view = local_view_host(sdf.weights());

                    auto field_view = local_view_host(local_field);
                    auto weights_view = local_view_host(local_weights);

                    auto node_mapping = mesh.node_mapping();

                    auto rr = out.range();
                    const int n_var = space.n_var();

                    for (SizeType i = 0; i < surface_field.size(); i++) {
                        const Scalar gg = surface_field_view.get(i);
                        // if (gg > cutoff) continue;

                        SizeType node = node_mapping[i];
                        // Offset for dof number
                        node *= n_var;

                        assert(node + 2 < local_is_contact.local_size());
                        assert(node + 2 < local_field.local_size());
                        assert(node + 2 < local_normals.local_size());

                        field_view.set(node, gg);
                        weights_view.set(node, surface_weights_view.get(i));
                    }
                }

                // As we add data to these vectors we have reset them to zero
                out.set(0);

                // Temp vector needs to be created
                space.create_vector(weights);
                weights.set(0);

                auto mode = utopia::ADD_MODE;
                // auto mode = utopia::OVERWRITE_MODE;
                space.local_to_global(local_field, out, mode);
                space.local_to_global(local_weights, weights, mode);

                // Divide by weight
                e_pseudo_inv(out, weights, device::Epsilon<Scalar>::value());
            } else if (!local_to_global.empty()) {
                // Scope for views
                {
                    auto surface_field_view = const_local_view_host(surface_field);
                    auto field_view = local_view_host(out);
                    auto node_mapping = mesh.node_mapping();

                    auto rr = out.range();
                    const SizeType local_size = out.local_size();
                    for (ptrdiff_t i = 0; i < mesh.n_local_nodes(); i++) {
                        const Scalar gg = surface_field_view.get(i);

                        SizeType node = node_mapping[i];
                        node = local_to_global(node, 0) - rr.begin();

                        if (node >= local_size || node < 0) {
                            // Skip ghost nodes!
                            continue;
                        }

                        assert(node + 2 < is_contact.local_size());
                        assert(node + 2 < sdf->data().local_size());
                        assert(node + 2 < normals->data().local_size());

                        field_view.set(node, gg);
                    }
                }

            } else {
                // Serial case
                auto surface_field_view = const_local_view_device(surface_field);
                auto field_view = local_view_device(out);

                auto node_mapping = mesh.node_mapping();

                const int n_var = space.n_var();
                for (ptrdiff_t i = 0; i < mesh.n_local_nodes(); i++) {
                    const Scalar gg = surface_field_view.get(i);
                    SizeType node = node_mapping[i];
                    // Offset for dof number
                    node *= n_var;

                    field_view.set(node, gg);
                }
            }
        }
    };

    template <class FunctionSpace>
    bool SDFResample<FunctionSpace>::apply(FunctionSpace &space, Vector &out) {
        impl_->resample_to_mesh_surface_of(space, out);
        return true;
    }

    template <class FunctionSpace>
    SDFResample<FunctionSpace>::SDFResample() : impl_(utopia::make_unique<Impl>()) {}

    template <class FunctionSpace>
    SDFResample<FunctionSpace>::~SDFResample() {}

    template <class FunctionSpace>
    void SDFResample<FunctionSpace>::read(Input &in) {
        impl_->read(in);
    }
}  // namespace utopia

#endif  // UTOPIA_SDF_RESAMPLE_IMPL_HPP

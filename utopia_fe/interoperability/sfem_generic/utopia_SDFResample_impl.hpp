#ifndef UTOPIA_SDF_RESAMPLE_IMPL_HPP
#define UTOPIA_SDF_RESAMPLE_IMPL_HPP

#include <petscvec.h>
#include "utopia_Epsilon.hpp"
#include "utopia_Input.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_sfem_Mesh.hpp"
#include "utopia_sfem_SDF.hpp"

#include "utopia_SDFResample.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_stk_DofMap.hpp"

#include <utopia_Enums.hpp>
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
        Scalar infinity{1e10};

        void read(Input &in) override {
            sdf.read(in);

            in.get("verbose", verbose);
            in.get("shift", shift);
            in.get("infinity", infinity);

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

            if (false)  //
            {
                static int export_mesh = 0;
                mesh.write("dbg_mesh" + std::to_string(export_mesh++));
            }

            Vector surface_field;
            sdf.to_mesh(mesh, surface_field);
            surface_field.shift(shift);

        
            if (verbose) {
                Scalar sum_gap = sum(surface_field);
                Scalar max_gap = max(surface_field);
                Scalar min_gap = min(surface_field);

                std::stringstream ss;
                ss << "Surface Mesh #nodes " << mesh.n_local_nodes() << "\n";
                ss << "sum_gap: " << sum_gap << "\n";
                ss << "min_gap: " << min_gap << "\n";
                ss << "max_gap: " << max_gap << "\n";
                space.comm().root_print(ss.str());
            }

            space.create_vector(out);

            auto &&local_to_global = space.dof_map().local_to_global();

            if (sdf.has_weights()) {
                Vector weights;
                space.create_vector(weights);

                weights.set(0);
                out.set(0);

                // Scope for views
                {
                    auto surface_field_view = const_local_view_host(surface_field);
                    auto surface_weights_view = local_view_host(sdf.weights());
                    auto node_mapping = mesh.node_mapping();

                    Write<Vector> w_out(out, utopia::GLOBAL_ADD);
                    Write<Vector> w_weights(weights, utopia::GLOBAL_ADD);

                    auto rr = out.range();
                    for (SizeType i = 0; i < surface_field.size(); i++) {
                        const Scalar value = surface_field_view.get(i);
                        const SizeType node = node_mapping[i];
                        const SizeType globalid = local_to_global(node, 0);

                        out.c_add(globalid, value);
                        weights.c_add(globalid, surface_weights_view.get(i));
                    }
                }

                // Divide by weight
                e_pseudo_inv(weights, weights, device::Epsilon<Scalar>::value());
                out = e_mul(out, weights);
            } else if (!local_to_global.empty()) {
                out.set(infinity);
                // Scope for views
                {
                    auto surface_field_view = const_local_view_host(surface_field);
                    auto node_mapping = mesh.node_mapping();
                    Write<Vector> w(out, utopia::GLOBAL_INSERT);

                    auto rr = out.range();
                    for (ptrdiff_t i = 0; i < node_mapping.size(); i++) {
                        const Scalar value = surface_field_view.get(i);

                        SizeType node = node_mapping[i];
                        SizeType globalid = local_to_global(node, 0);
                        out.c_set(globalid, value);
                    }
                }

            } else {
                out.set(infinity);
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

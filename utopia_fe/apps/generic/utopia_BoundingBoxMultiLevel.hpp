#ifndef UTOPIA_BOUNDING_BOX_MULTILEVEL_HPP
#define UTOPIA_BOUNDING_BOX_MULTILEVEL_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

namespace utopia {

    template <class Mesh>
    class BoundingBoxMultiLevelMeshGenerator : public Configurable {
    public:
        using AABB = typename Traits<Mesh>::AABB;
        using Scalar = typename Traits<Mesh>::Scalar;

        void read(Input &in) override {
            in.get("verbose", verbose_);
            in.get("elem_type", elem_type_);
        }

        bool generate_meshes(Mesh &fine_mesh, const int n_levels, std::vector<std::shared_ptr<Mesh>> &meshes) {
            assert(n_levels > 0);
            if (n_levels == 0) return false;

            AABB box;
            fine_mesh.bounding_box(box);

            typename AABB::Point r = box.max;

            const int spatial_dim = fine_mesh.spatial_dimension();

            Scalar max_r = 0.0;
            for (int d = 0; d < spatial_dim; ++d) {
                r[d] -= box.min[d];
                max_r = std::max(max_r, r[d]);
            }

            SizeType n_elements = fine_mesh.n_elements();
            SizeType n_nodes = fine_mesh.n_nodes();
            SizeType min_segments = 2;
            SizeType computed_segments =
                std::pow(Scalar(n_nodes) / std::pow(2, spatial_dim * n_levels), Scalar(1. / spatial_dim));

            const int n_segments = std::max(min_segments, computed_segments);

            std::array<SizeType, 3> n = {0, 0, 0};

            for (int d = 0; d < spatial_dim; ++d) {
                double aspect_ratio = r[d] / max_r;
                n[d] = n_segments * aspect_ratio;
            }

            if (verbose_) {
                std::stringstream ss;
                ss << "Generating hierarchy from mesh with " << n_nodes << " nodes and " << n_elements
                   << " elements.\n";

                fine_mesh.comm().root_print(ss.str());
            }

            for (int l = 0; l < n_levels; ++l) {
                auto mesh = std::make_shared<Mesh>(fine_mesh.comm());
                mesh->box(box, n[0], n[1], n[2], elem_type_);
                meshes.push_back(mesh);

                if (verbose_) {
                    std::stringstream ss;
                    ss << "Generated mesh at level " << l << " with " << n[0] << " x " << n[1] << " x " << n[2]
                       << " intervals -> " << mesh->n_elements() << " elements\n";

                    fine_mesh.comm().root_print(ss.str());
                }

                for (int d = 0; d < spatial_dim; ++d) {
                    n[d] *= 2;
                }
            }

            return true;
        }

    private:
        bool verbose_{true};
        std::string elem_type_{"HEX8"};
    };

    template <class FunctionSpace>
    class BoundingBoxMultiLevelFunctionSpaceGenerator : public Configurable {
    public:
        using Mesh = typename Traits<FunctionSpace>::Mesh;
        using AABB = typename Traits<Mesh>::AABB;
        using Scalar = typename Traits<FunctionSpace>::Scalar;

        void read(Input &in) override { mesh_generator_.read(in); }

        bool generate_spaces(FunctionSpace &fine_space,
                             const int n_levels,
                             std::vector<std::shared_ptr<FunctionSpace>> &spaces) {
            std::vector<std::shared_ptr<Mesh>> meshes;
            if (!mesh_generator_.generate_meshes(fine_space.mesh(), n_levels, meshes)) {
                return false;
            }

            spaces.clear();

            for (auto &m_ptr : meshes) {
                auto space = std::make_shared<FunctionSpace>(m_ptr);
                space->copy_meta_info_from(fine_space);
                space->initialize();
                spaces.push_back(space);
            }

            return true;
        }

        BoundingBoxMultiLevelMeshGenerator<Mesh> mesh_generator_;
    };

}  // namespace utopia

#endif  // UTOPIA_BOUNDING_BOX_MULTILEVEL_HPP
#include "utopia_libmesh_Obstacle.hpp"

#include "moonolith_obstacle.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_make_unique.hpp"

#include "moonolith_affine_transform.hpp"
#include "moonolith_assign_functions.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_elem_quad.hpp"
#include "moonolith_elem_segment.hpp"
#include "moonolith_elem_shape.hpp"
#include "moonolith_elem_triangle.hpp"
#include "moonolith_matlab_scripter.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_sparse_matrix.hpp"

// #include "utopia_ConvertMesh.hpp"

#include "utopia_LibMeshToMoonolithConvertions.hpp"
#include "utopia_TransferUtils.hpp"

namespace utopia {

    namespace libmesh {

        class Obstacle::Impl : public Configurable {
        public:
            using Mesh_t = utopia::moonolith::Mesh;
            using FunctionSpace_t = utopia::moonolith::FunctionSpace;
            using Obstacle_t = utopia::moonolith::Obstacle;

            bool assemble(const FunctionSpace &in_space) {
                switch (in_space.mesh().spatial_dimension()) {
                    case 2: {
                        assert(space.raw_type<2>());
                        extract_trace_space(in_space.mesh().raw_type(),
                                            in_space.raw_type_dof_map(),
                                            obstacle.params().variable_number,
                                            *space.raw_type<2>(),
                                            {});
                        break;
                    }
                    case 3: {
                        assert(space.raw_type<3>());
                        extract_trace_space(in_space.mesh().raw_type(),
                                            in_space.raw_type_dof_map(),
                                            obstacle.params().variable_number,
                                            *space.raw_type<3>(),
                                            {});
                        break;
                    }
                    default: {
                        assert(false);
                        return false;
                    }
                }

                return obstacle.assemble(space);
            }

            bool init_obstacle(const Mesh &mesh) {
                if (mesh.manifold_dimension() == mesh.spatial_dimension() - 1) {
                    assert(false && "IMPLEMENT ME");
                    // switch (mesh.spatial_dimension()) {
                    //     case 2: {
                    //         assert(obstacle_mesh.raw_type<2>());
                    //         convert_mesh(mesh.raw_type(), *obstacle_mesh.raw_type<2>(), {});
                    //         break;
                    //     }
                    //     case 3: {
                    //         assert(obstacle_mesh.raw_type<3>());
                    //         convert_mesh(mesh.raw_type(), *obstacle_mesh.raw_type<3>(), {});
                    //         break;
                    //     }
                    //     default: {
                    //         assert(false);
                    //         return false;
                    //     }
                    // }

                } else {
                    assert(mesh.manifold_dimension() == mesh.spatial_dimension() - 1);
                    // FIXME THIS SHOULD GO IN SOME FUNCTION
                    switch (mesh.spatial_dimension()) {
                        case 2: {
                            assert(obstacle_mesh.raw_type<2>());
                            extract_surface<2>(mesh.raw_type(), *obstacle_mesh.raw_type<2>(), {});
                            break;
                        }
                        case 3: {
                            assert(obstacle_mesh.raw_type<3>());
                            extract_surface<3>(mesh.raw_type(), *obstacle_mesh.raw_type<3>(), {});
                            break;
                        }
                        default: {
                            assert(false);
                            return false;
                        }
                    }
                }

                return obstacle.init_obstacle(obstacle_mesh);
            }

            void read(Input &in) override { obstacle.read(in); }

            FunctionSpace_t space;
            Mesh_t obstacle_mesh;
            Obstacle_t obstacle;
        };

        Obstacle::Obstacle() : impl_(utopia::make_unique<Impl>()) {}
        Obstacle::~Obstacle() {}

        bool Obstacle::assemble(const FunctionSpace &space) { return impl_->assemble(space); }

        bool Obstacle::init_obstacle(const Mesh &obstacle_mesh) { return impl_->init_obstacle(obstacle_mesh); }

        const Obstacle::Vector &Obstacle::gap() const { return impl_->obstacle.gap(); }
        const Obstacle::Vector &Obstacle::is_contact() const { return impl_->obstacle.is_contact(); }

        void Obstacle::set_params(const Params &params) { impl_->obstacle.set_params(params); }

        void Obstacle::transform(const Matrix &in, Matrix &out) { impl_->obstacle.transform(in, out); }

        void Obstacle::transform(const Vector &in, Vector &out) { impl_->obstacle.transform(in, out); }

        void Obstacle::inverse_transform(const Vector &in, Vector &out) { impl_->obstacle.inverse_transform(in, out); }

    }  // namespace libmesh
}  // namespace utopia

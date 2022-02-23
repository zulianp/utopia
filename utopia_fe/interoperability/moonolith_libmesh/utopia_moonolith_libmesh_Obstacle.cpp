#include "utopia_moonolith_libmesh_Obstacle.hpp"

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

#include "utopia_TransferUtils.hpp"
#include "utopia_moonolith_libmesh_Convert.hpp"

namespace utopia {

    namespace libmesh {

        class Obstacle::Impl : public Configurable {
        public:
            using Mesh_t = utopia::moonolith::Mesh;
            using FunctionSpace_t = utopia::moonolith::FunctionSpace;
            using Obstacle_t = utopia::moonolith::Obstacle;

            bool assemble(const FunctionSpace &in_space) {
                // FIXME THIS SHOULD GO IN SOME REUSABLE FUNCTION

                switch (in_space.mesh().spatial_dimension()) {
                    case 2: {
                        using Mesh2_t = ::moonolith::Mesh<Scalar, 2>;
                        using Space2_t = ::moonolith::FunctionSpace<Mesh2_t>;

                        auto mesh = std::make_shared<Mesh2_t>(in_space.comm().raw_comm());
                        auto temp_space = std::make_shared<Space2_t>(mesh);

                        extract_trace_space(in_space.mesh().raw_type(),
                                            in_space.raw_type_dof_map(),
                                            obstacle.params().variable_number,
                                            *temp_space,
                                            {});

                        space.wrap(temp_space);

                        break;
                    }
                    case 3: {
                        using Mesh3_t = ::moonolith::Mesh<Scalar, 3>;
                        using Space3_t = ::moonolith::FunctionSpace<Mesh3_t>;

                        auto mesh = std::make_shared<Mesh3_t>(in_space.comm().raw_comm());
                        auto temp_space = std::make_shared<Space3_t>(mesh);

                        extract_trace_space(in_space.mesh().raw_type(),
                                            in_space.raw_type_dof_map(),
                                            obstacle.params().variable_number,
                                            *temp_space,
                                            {});

                        space.wrap(temp_space);
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
                    assert(mesh.manifold_dimension() == mesh.spatial_dimension());
                    // FIXME THIS SHOULD GO IN SOME REUSABLE FUNCTION
                    switch (mesh.spatial_dimension()) {
                        case 2: {
                            auto m_mesh = std::make_shared<::moonolith::Mesh<Scalar, 2>>(mesh.comm().raw_comm());
                            extract_surface<2>(mesh.raw_type(), *m_mesh, {});
                            obstacle_mesh.wrap(m_mesh);
                            break;
                        }
                        case 3: {
                            auto m_mesh = std::make_shared<::moonolith::Mesh<Scalar, 3>>(mesh.comm().raw_comm());
                            extract_surface<3>(mesh.raw_type(), *m_mesh, {});
                            obstacle_mesh.wrap(m_mesh);
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

        void Obstacle::read(Input &in) { impl_->obstacle.read(in); }
        void Obstacle::describe(std::ostream &os) const { impl_->obstacle.describe(os); }

        bool Obstacle::assemble(FunctionSpace &space) { return impl_->assemble(space); }

        bool Obstacle::init_obstacle(const Mesh &obstacle_mesh) { return impl_->init_obstacle(obstacle_mesh); }

        const Obstacle::Vector &Obstacle::gap() const { return impl_->obstacle.gap(); }
        const Obstacle::Vector &Obstacle::is_contact() const { return impl_->obstacle.is_contact(); }
        const Obstacle::Vector &Obstacle::normals() const { return impl_->obstacle.normals(); }

        void Obstacle::set_params(const Params &params) { impl_->obstacle.set_params(params); }

        void Obstacle::transform(const Matrix &in, Matrix &out) { impl_->obstacle.transform(in, out); }

        void Obstacle::transform(const Vector &in, Vector &out) { impl_->obstacle.transform(in, out); }

        void Obstacle::inverse_transform(const Vector &in, Vector &out) { impl_->obstacle.inverse_transform(in, out); }

    }  // namespace libmesh
}  // namespace utopia

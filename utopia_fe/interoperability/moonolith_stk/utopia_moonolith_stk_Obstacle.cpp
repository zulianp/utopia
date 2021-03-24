#include "utopia_moonolith_stk_Obstacle.hpp"

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

namespace utopia {

    namespace stk {

        class Obstacle::Impl : public Configurable {
        public:
            using Mesh_t = utopia::moonolith::Mesh;
            using FunctionSpace_t = utopia::moonolith::FunctionSpace;
            using Obstacle_t = utopia::moonolith::Obstacle;

            bool assemble(const FunctionSpace &in_space) {}

            bool init_obstacle(const Mesh &mesh) {}

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
        const Obstacle::Vector &Obstacle::normals() const { return impl_->obstacle.normals(); }

        void Obstacle::set_params(const Params &params) { impl_->obstacle.set_params(params); }

        void Obstacle::transform(const Matrix &in, Matrix &out) { impl_->obstacle.transform(in, out); }

        void Obstacle::transform(const Vector &in, Vector &out) { impl_->obstacle.transform(in, out); }

        void Obstacle::inverse_transform(const Vector &in, Vector &out) { impl_->obstacle.inverse_transform(in, out); }

    }  // namespace stk
}  // namespace utopia

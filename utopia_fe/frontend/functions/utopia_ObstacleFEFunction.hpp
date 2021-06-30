#ifndef UTOPIA_OBSTACLEF_EF_UNCTION_HPP
#define UTOPIA_OBSTACLEF_EF_UNCTION_HPP

#include "utopia_BoxConstrainedFEFunction.hpp"
#include "utopia_fe_Core.hpp"

#include "utopia_ImplicitObstacle_impl.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ObstacleFEFunction final : public BoxConstrainedFEFunction<FunctionSpace> {
    public:
        using Super = utopia::BoxConstrainedFEFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;

        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;

        // Use specialized components for function space
        using Obstacle_t = utopia::Obstacle<FunctionSpace>;
        using ImplicitObstacle_t = utopia::ImplicitObstacle<FunctionSpace>;

        bool has_nonlinear_constraints() const override { return !linear_obstacle_; }

        bool update_constraints(const Vector_t &x) {
            utopia::out() << "update_constraints\n";

            this->space()->displace(x);
            bool ok = obstacle_->assemble(*this->space());

            if (debug_) {
                static int iter_debug = 0;
                ouput_debug_data(iter_debug++, x);
            }

            this->space()->displace(-x);
            return ok;
        }

        bool has_transformation() const override { return true; }

        void transform(const Vector_t &x, Vector_t &x_constrained) override { obstacle_->transform(x, x_constrained); }

        void inverse_transform(const Vector_t &x_constrained, Vector_t &x) override {
            obstacle_->inverse_transform(x_constrained, x);
        }

        void transform(const Matrix_t &H, Matrix_t &H_constrained) override { obstacle_->transform(H, H_constrained); }

        bool constraints_gradient(const Vector_t &x, BoxConstraints<Vector_t> &box) override {
            update_constraints(x);

            if (!box.upper_bound()) {
                box.upper_bound() = std::make_shared<Vector_t>();
            }

            (*box.upper_bound()) = obstacle_->gap();
            return true;
        }

        ObstacleFEFunction(const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &unconstrained)
            : Super(unconstrained) {}

        void ouput_debug_data(const Size_t iter_debug, const Vector_t &) const {
            // Output extras
            {
                Vector_t gap_zeroed = e_mul(obstacle_->is_contact(), obstacle_->gap());
                Vector_t rays = obstacle_->normals();

                {
                    auto rays_view = local_view_device(rays);
                    auto gap_view = const_local_view_device(obstacle_->gap());

                    int dim = this->space()->mesh().spatial_dimension();

                    parallel_for(
                        local_range_device(rays), UTOPIA_LAMBDA(const Size_t i) {
                            Size_t block = i / dim;
                            Size_t g_idx = dim * block;

                            auto ri = rays_view.get(i);
                            ri *= gap_view.get(g_idx);
                            rays_view.set(i, ri);
                        });
                }

                this->space()->write("rays_" + std::to_string(iter_debug) + ".e", rays);
            }
        }

        void read(Input &in) override {
            Super::read(in);

            in.get("linear_obstacle", linear_obstacle_);
            in.get("debug", debug_);

            if (!obstacle_) {
                auto implicit_obstacle = std::make_shared<ImplicitObstacle_t>();
                in.get("implicit_obstacle", *implicit_obstacle);

                obstacle_ = std::make_shared<Obstacle_t>();
                in.require("obstacle", params_);
                // Must be created for every process independently and the same
                Mesh_t obstacle_mesh(Communicator_t::self());
                in.require("obstacle", obstacle_mesh);

                obstacle_->set_params(params_);
                obstacle_->init_obstacle(obstacle_mesh);

                bool export_obstacle = false;
                in.get("export_obstacle", export_obstacle);

                if (export_obstacle) {
                    if (this->space()->comm().rank() == 0) {
                        obstacle_mesh.write("obstacle.e");
                    }
                }
            }
        }

        // bool report_solution(const Vector_t &x) override {
        //     if (!Super::report_solution(x)) {
        //         return false;
        //     }

        //     if (debug_) {
        //         static int iter_debug = 0;
        //         ouput_debug_data(iter_debug++, x);
        //     }

        //     return true;
        // }

    private:
        typename Obstacle_t::Params params_;
        std::shared_ptr<Obstacle_t> obstacle_;
        bool linear_obstacle_{false};
        bool debug_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_FE_FUNCTION_HPP

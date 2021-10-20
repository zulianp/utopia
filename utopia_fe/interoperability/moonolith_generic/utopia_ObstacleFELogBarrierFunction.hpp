#ifndef UTOPIA_OBSTACLE_FE_LOG_BARRIER_FUNCTION_HPP
#define UTOPIA_OBSTACLE_FE_LOG_BARRIER_FUNCTION_HPP

#include "utopia_LogBarrierFunctionFactory.hpp"

#include "utopia_BoxConstrainedFEFunction.hpp"
#include "utopia_IObstacle.hpp"
#include "utopia_fe_Core.hpp"

#include "utopia_AnalyticObstacle_impl.hpp"
#include "utopia_ImplicitObstacle_impl.hpp"

#include "utopia_ObstacleFactory.hpp"

#include "utopia_IsNotSupported.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ObstacleFELogBarrierFunction final : public FEFunctionInterface<FunctionSpace> {
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
        using AnalyticObstacle_t = utopia::AnalyticObstacle<FunctionSpace>;

        // bool has_nonlinear_constraints() const override { return !linear_obstacle_; }

        // inline std::shared_ptr<Vector_t> selection() override {
        //     auto s = std::make_shared<Vector_t>(obstacle_->is_contact());
        //     // this->space()->apply_zero_constraints(*s);
        //     return s;
        // }

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

        ObstacleFELogBarrierFunction(const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &unconstrained)
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

            in.get("debug", debug_);

            if (!obstacle_) {
                std::string type;
                in.get("obstacle",
                       [&](Input &node) { obstacle_ = ObstacleFactory<FunctionSpace>::new_obstacle(node); });
            }

            std::string function_type;
            function_ = LogBarrierFunctionFactory<Matrix, Vector>::new_log_barrier_function(function_type);
            function_->read(in);
        }

    private:
        std::shared_ptr<IObstacle<FunctionSpace>> obstacle_;
        std::shared_ptr<LogBarrierFunctionBase> function_;
        bool debug_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_FE_LOG_BARRIER_FUNCTION_HPP
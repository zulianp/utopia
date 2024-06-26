#ifndef UTOPIA_OBSTACLE_NEWMARK_HPP
#define UTOPIA_OBSTACLE_NEWMARK_HPP

#include "utopia_LogBarrierFactory.hpp"

#include "utopia_BoxConstrainedFEFunction.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_ContactInterface.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_fe_Core.hpp"

#include "utopia_AnalyticObstacle_impl.hpp"
#include "utopia_ImplicitObstacle_impl.hpp"

#include "utopia_ContactFactory.hpp"

#include "utopia_IsNotSupported.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ObstacleNewmark final : public NewmarkIntegrator<FunctionSpace>,
                                  public LSStrategy<typename Traits<FunctionSpace>::Vector> {
    public:
        using Super = utopia::NewmarkIntegrator<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Layout_t = typename Traits<FunctionSpace>::Layout;

        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;
        using Penalty = utopia::Penalty<Matrix_t, Vector_t>;

        ObstacleNewmark(const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &unconstrained)
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
            in.get("debug_from_iteration", debug_from_iteration_);
            in.get("trivial_obstacle", trivial_obstacle_);
            in.get("allow_projection", allow_projection_);
            in.get("enable_line_search", enable_line_search_);

            if (!obstacle_) {
                std::string type;
                in.get("obstacle", [&](Input &node) { obstacle_ = ContactFactory<FunctionSpace>::new_obstacle(node); });
            }

            if (!obstacle_) {
                std::string type;
                in.get("contact", [&](Input &node) { obstacle_ = ContactFactory<FunctionSpace>::new_contact(node); });
            }

            ////////////////////////////////////////////////////////////////////////////////
            std::string function_type;
            in.get("function_type", function_type);
            barrier_ = PenaltyFactory<Matrix_t, Vector_t>::new_penalty(function_type);
            barrier_->read(in);
        }

        bool setup_IVP(Vector_t &x) override {
            update_constraints(x);
            barrier_->reset();

            return Super::setup_IVP(x);
        }

        bool time_derivative(const Vector_t &x, Vector_t &dfdt) const override {
            return Super::time_derivative(x, dfdt);
        }

        bool update_IVP(const Vector_t &x) override {
            if (!trivial_obstacle_) {
                update_constraints(x);
            }

            if (x.has_nan_or_inf()) {
                this->~ObstacleNewmark();
                x.comm().barrier();
                assert(false);
                Utopia::Abort("update_IVP: NaN detected!");
            }

            barrier_->reset();
            return Super::update_IVP(x);
        }

        void barrier_hessian(const Vector_t &x, Matrix_t &H) const {
            if (barrier_) {
                Vector_t barrier_temp(layout(x), 0);
                Matrix_t barrier_H, temp_H;

                if (trivial_obstacle_) {
                    obstacle_->transform(x, barrier_temp);
                } else {
                    Vector_t delta_x = x - this->x_old();
                    obstacle_->transform(delta_x, barrier_temp);
                }

                barrier_H.identity(square_matrix_layout(layout(x)), 0);

                barrier_->hessian(barrier_temp, barrier_H);

                obstacle_->transform(barrier_H, temp_H);

                // Remove contribution from boundary conditions
                this->space()->apply_constraints(temp_H, 0);
                H += temp_H;
            }
        }

        void barrier_gradient(const Vector_t &x, Vector_t &g) const {
            if (barrier_) {
                Vector_t barrier_temp(layout(x), 0);
                Vector_t barrier_g(layout(g), 0);

                if (trivial_obstacle_) {
                    obstacle_->transform(x, barrier_temp);
                } else {
                    // Use barrier_g as a temporary for the delta_x
                    barrier_g = x - this->x_old();
                    obstacle_->transform(barrier_g, barrier_temp);

                    // Reset barrier_g to zero
                    barrier_g.set(0.);
                }

                barrier_->gradient(barrier_temp, barrier_g);

                obstacle_->inverse_transform(barrier_g, barrier_temp);
                g += barrier_temp;
                // g -= barrier_temp;

                if (g.has_nan_or_inf()) {
                    this->~ObstacleNewmark();
                    g.comm().barrier();
                    assert(false);
                    Utopia::Abort("barrier_gradient: NaN detected!");
                }
            }
        }

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            Super::integrate_gradient(x, g);
            barrier_gradient(x, g);
        }

        void integrate_hessian(const Vector_t &x, Matrix_t &H) const override {
            UTOPIA_TRACE_REGION_BEGIN("ObstacleNewmark::integrate_hessian");
            Super::integrate_hessian(x, H);
            barrier_hessian(x, H);
            UTOPIA_TRACE_REGION_END("ObstacleNewmark::integrate_hessian");
        }

        ////////////////////////////////////////////////////////////////////////////////

        const Vector_t &solution() const override { return this->x_old(); }
        bool report_solution(const Vector_t &x) override { return Super::report_solution(x); }

        // bool update(const Vector_t &x) override {
        //     if (barrier_) {
        //         barrier_->update(x);
        //     }

        //     return true;
        // }

        bool update(const Vector_t &x) override {
            if (barrier_) {
                Vector_t x_obs;

                if (trivial_obstacle_) {
                    obstacle_->transform(x, x_obs);
                } else {
                    Vector_t x_temp = x - this->x_old();
                    obstacle_->transform(x_temp, x_obs);
                }

                barrier_->update(x_obs);
            }

            return true;
        }

        bool project_onto_feasibile_region(Vector_t &x) const final {
            bool ok = true;

            if (!allow_projection_) return ok;

            if (barrier_) {
                Vector_t barrier_x(layout(x), 0);
                Vector_t delta_x;

                if (trivial_obstacle_) {
                    obstacle_->transform(x, barrier_x);
                } else {
                    delta_x = x - this->x_old();
                    obstacle_->transform(delta_x, barrier_x);
                }

                ok = barrier_->project_onto_feasibile_region(barrier_x);

                if (trivial_obstacle_) {
                    obstacle_->inverse_transform(barrier_x, x);
                } else {
                    obstacle_->inverse_transform(barrier_x, delta_x);
                    x = this->x_old() + delta_x;
                }
            }

            return ok;
        }

        bool get_alpha(FunctionBase<Vector_t> &fun,
                       const Vector_t &g,
                       const Vector_t &x,
                       const Vector_t &correction,
                       Scalar_t &alpha) override {
            if (line_search_) {
                if (trivial_obstacle_) {
                    Vector_t zero(layout(x), 0.);
                    return line_search_->get_alpha(fun, g, zero, correction, alpha);
                } else {
                    return line_search_->get_alpha(fun, g, x, correction, alpha);
                }
            } else {
                alpha = 0.99;
                return false;
            }
        }

        bool get_alpha(LeastSquaresFunctionBase<Vector_t> &fun,
                       const Vector_t &g,
                       const Vector_t &x,
                       const Vector_t &correction,
                       Scalar_t &alpha) override {
            if (line_search_) {
                if (trivial_obstacle_) {
                    Vector_t zero(layout(x), 0.);
                    return line_search_->get_alpha(fun, g, zero, correction, alpha);
                } else {
                    return line_search_->get_alpha(fun, g, x, correction, alpha);
                }

            } else {
                alpha = 0.99;
                return false;
            }
        }

        void init_memory(const Layout_t & /*layout*/) override {}

        inline std::shared_ptr<LSStrategy<Vector_t>> line_search() override {
            if (line_search_) {
                return make_ref(*this);
            } else {
                return nullptr;
            }
        }

    private:
        std::shared_ptr<ContactInterface<FunctionSpace>> obstacle_;
        std::shared_ptr<Penalty> barrier_;
        bool debug_{false};
        int debug_from_iteration_{0};
        bool trivial_obstacle_{false};
        bool allow_projection_{true};

        bool enable_line_search_{false};
        std::shared_ptr<LineSearchBoxProjection<Vector_t>> line_search_;

        bool update_constraints(const Vector_t &x) {
            utopia::out() << "update_constraints\n";

            this->space()->displace(x);
            bool ok = obstacle_->assemble(*this->space());

            if (debug_) {
                static int iter_debug = 0;
                if (iter_debug >= debug_from_iteration_) {
                    ouput_debug_data(iter_debug, x);
                }

                iter_debug++;
            }

            this->space()->displace(-x);

            if (!barrier_) {
                x.comm().root_print("[Warning] no barrier!");
                barrier_ = std::make_shared<LogBarrier<Matrix_t, Vector_t>>();
            }

            auto box =
                std::make_shared<BoxConstraints<Vector_t>>(nullptr, std::make_shared<Vector_t>(obstacle_->gap()));
            barrier_->set_box_constraints(box);

            barrier_->set_selection(std::make_shared<Vector_t>(obstacle_->is_contact()));
            barrier_->set_scaling_matrix(obstacle_->mass_matrix());

            if (enable_line_search_) {
                if (!line_search_) {
                    line_search_ = std::make_shared<LineSearchBoxProjection<Vector_t>>(box, make_ref(this->x_old()));
                } else {
                    line_search_->set_box_constraints(box);
                    line_search_->set_offset_vector(make_ref(this->x_old()));
                }

                line_search_->set_dumping(0.98);

                auto trafo = obstacle_->orthogonal_transformation();
                assert(trafo);
                if (!trafo) {
                    Utopia::Abort(
                        "ObstacleVelocityNewmark:update_constraints: orthogonal_transformation is mandatory for "
                        "line_search!");
                }

                line_search_->set_transform(trafo);
            }

            return ok;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_NEWMARK_HPP

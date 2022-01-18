#ifndef UTOPIA_OBSTACLE_STABILIZED_VELOCITY_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_OBSTACLE_STABILIZED_VELOCITY_NEWMARK_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"
#include "utopia_IObstacle.hpp"
#include "utopia_ObstacleFactory.hpp"

#include "utopia_LineSearchBoxProjection.hpp"
#include "utopia_LogBarrierFactory.hpp"
#include "utopia_ObstacleStabilizedNewmark.hpp"

#include <utility>

namespace utopia {

    // https://en.wikipedia.org/wiki/Newmark-beta_method
    // Unconditionally Stable gamma = 0.5, beta = 0.25
    template <class FunctionSpace>
    class ObstacleStabilizedVelocityNewmark : public TimeDependentFunction<FunctionSpace>,
                                              public ObstacleDependentFunction<FunctionSpace>,
                                              public LSStrategy<typename Traits<FunctionSpace>::Vector> {
    public:
        using Super = utopia::TimeDependentFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Layout_t = typename Traits<FunctionSpace>::Layout;

        using LogBarrierBase = utopia::LogBarrierBase<Matrix_t, Vector_t>;

        void read(Input &in) override {
            Super::read(in);

            in.get("debug", debug_);
            in.get("debug_from_iteration", debug_from_iteration_);
            in.get("stabilized_formulation", stabilized_formulation_);
            in.get("trivial_obstacle", trivial_obstacle_);
            in.get("zero_initial_guess", zero_initial_guess_);
            in.get("enable_line_search", enable_line_search_);

            if (!obstacle_) {
                std::string type;
                in.get("obstacle",
                       [&](Input &node) { obstacle_ = ObstacleFactory<FunctionSpace>::new_obstacle(node); });
            }

            if (stabilized_formulation_) {
                utopia::out() << "Using stabilized formulation!\n";
            }

            ////////////////////////////////////////////////////////////////////////////////
            std::string function_type;
            in.get("function_type", function_type);
            barrier_ = LogBarrierFactory<Matrix_t, Vector_t>::new_log_barrier(function_type);
            barrier_->read(in);

            bool use_barrier_mass_scaling = false;
            in.get("use_barrier_mass_scaling", use_barrier_mass_scaling);

            if (use_barrier_mass_scaling) {
                barrier_->set_scaling_matrix(obstacle_->mass_matrix());
            }
        }

        bool update_constraints(const Vector_t &x) {
            utopia::out() << "ObstacleStabilizedVelocityNewmark::update_constraints\n";

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

            if (enable_line_search_) {
                if (!line_search_) {
                    line_search_ = std::make_shared<LineSearchBoxProjection<Vector_t>>(box, make_ref(this->x_old()));
                } else {
                    line_search_->set_box_constraints(box);
                    line_search_->set_offset_vector(make_ref(this->x_old()));
                }

                auto trafo = obstacle_->orthogonal_transformation();
                assert(trafo);
                if (!trafo) {
                    Utopia::Abort(
                        "ObstacleStabilizedVelocityNewmark:update_constraints: orthogonal_transformation is mandatory "
                        "for "
                        "line_search!");
                }

                line_search_->set_transform(trafo);
            }

            return ok;
        }

        bool get_alpha(FunctionBase<Vector_t> &fun,
                       const Vector_t &g,
                       const Vector_t &velocity,
                       const Vector_t &correction,
                       Scalar_t &alpha) override {
            if (line_search_) {
                Vector_t c;
                update_x(velocity, c);
                c -= this->x_old();

                if (trivial_obstacle_) {
                    Vector_t zero(layout(c), 0.);
                    return line_search_->get_alpha(fun, g, zero, correction, alpha);
                } else {
                    return line_search_->get_alpha(fun, g, this->x_old(), correction, alpha);
                }
            } else {
                alpha = 1.;
                return false;
            }
        }

        bool get_alpha(LeastSquaresFunctionBase<Vector_t> &fun,
                       const Vector_t &g,
                       const Vector_t &velocity,
                       const Vector_t &correction,
                       Scalar_t &alpha) override {
            if (line_search_) {
                Vector_t c;
                update_x(velocity, c);
                c -= this->x_old();

                if (trivial_obstacle_) {
                    Vector_t zero(layout(c), 0.);
                    return line_search_->get_alpha(fun, g, zero, correction, alpha);
                } else {
                    return line_search_->get_alpha(fun, g, this->x_old(), correction, alpha);
                }

            } else {
                alpha = 1.;
                return false;
            }
        }

        void init_memory(const Layout_t & /*layout*/) override {}

        void update_predictor() {
            assert(obstacle_);

            if (stabilized_formulation_) {
                const Scalar_t dt = this->delta_time();
                Vector_t temp;
                obstacle_->transform(velocity_old_, temp);
                temp *= dt;
                temp = utopia::min(temp, obstacle_->gap());

                obstacle_->inverse_transform(temp, predictor_);
            } else {
                predictor_ = this->delta_time() * velocity_old_;
            }

            predictor_ += x_old_;
        }

        bool setup_IVP(Vector_t &x) override {
            if (!this->assemble_mass_matrix()) {
                return false;
            }

            assert(this->mass_matrix());
            Scalar_t sum_mm = sum(*this->mass_matrix());
            has_zero_density_ = sum_mm == 0.0;

            auto vlo = layout(x);

            x_old_ = x;
            velocity_old_.zeros(vlo);
            force_old_.zeros(vlo);

            update_constraints(x);

            update_predictor();
            return true;
        }

        bool update_IVP(const Vector_t &velocity) override {
            Vector_t x;
            update_x(velocity, x);

            if (x.has_nan_or_inf()) {
                this->~ObstacleStabilizedVelocityNewmark();
                x.comm().barrier();
                assert(false);
                Utopia::Abort("update_IVP: NaN detected!");
            }

            this->function()->gradient(x, force_old_);

            // Store current solution
            x_old_ = x;
            velocity_old_ = velocity;

            if (!trivial_obstacle_) {
                update_constraints(x);
            }

            bool ok = Super::update_IVP(x);
            update_predictor();
            barrier_->reset();
            return ok;
        }

        bool time_derivative(const Vector_t &x, Vector_t &velocity) const override {
            const Scalar_t dt = this->delta_time();
            velocity = velocity_old_;
            velocity += (2 / dt) * (x - predictor_);
            return true;
        }

        template <class... Args>
        ObstacleStabilizedVelocityNewmark(Args &&... args) : Super(std::forward<Args>(args)...) {}

        virtual ~ObstacleStabilizedVelocityNewmark() = default;

        bool gradient(const Vector_t &velocity, Vector_t &g) const override {
            Vector_t x;
            update_x(velocity, x);
            if (!this->function()->gradient(x, g)) {
                return false;
            }

            integrate_gradient(x, g);
            return true;
        }

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();

            if (!has_zero_density_) {
                Vector_t mom = (x - predictor_);
                mom *= (4.0 / dt2);
                g += (*this->mass_matrix()) * mom;
                g += force_old_;
            }

            barrier_gradient(x, g);

            if (this->must_apply_constraints_to_assembled()) {
                this->space()->apply_zero_constraints(g);
            }
        }

        bool hessian(const Vector_t &velocity, Matrix_t &H) const override {
            Vector_t x;
            update_x(velocity, x);
            return Super::hessian(x, H);
        }

        void integrate_hessian(const Vector_t &x, Matrix_t &H) const override {
            if (!has_zero_density_) {
                const Scalar_t dt2 = this->delta_time() * this->delta_time();
                H += (4. / dt2) * (*this->mass_matrix());
            }

            barrier_hessian(x, H);
            H *= (this->delta_time() / 2);

            this->space()->apply_constraints(H);
        }

        bool hessian_and_gradient(const Vector_t &velocity, Matrix_t &H, Vector_t &g) const override {
            Vector_t x;
            update_x(velocity, x);
            return Super::hessian_and_gradient(x, H, g);
        }

        bool update(const Vector_t & /*x*/) override {
            if (barrier_) {
                barrier_->update_barrier();
            }

            return true;
        }

        inline Vector_t &x_old() { return x_old_; }
        inline const Vector_t &x_old() const { return x_old_; }

        const Vector_t &solution() const override { return x_old(); }
        const Vector_t &velocity() const { return velocity_old_; }

        bool report_solution(const Vector_t &) override { return Super::report_solution(solution()); }

        bool set_initial_condition(const Vector_t &x) override {
            x_old_ = x;
            return true;
        }

        inline void set_obstacle(const std::shared_ptr<IObstacle<FunctionSpace>> obstacle) override {
            obstacle_ = obstacle;
        }

        bool project_onto_feasibile_region(Vector_t &velocity) const final {
            bool ok = true;
            if (barrier_) {
                Vector_t x;
                update_x(velocity, x);

                project_x_onto_feasibile_region(x);

                this->time_derivative(x, velocity);
            }

            return ok;
        }

        bool project_x_onto_feasibile_region(Vector_t &x) const {
            bool ok = true;

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

        void initial_guess_for_solver(Vector_t &velocity) override {
            if (zero_initial_guess_) {
                velocity.set(0.);
                return;
            }

            Vector_t x = this->x_old();
            update_x(velocity, x);
            x -= this->x_old();

            Scalar_t alpha = 0.99;
            // if (line_search_) {
            //     alpha = line_search_->compute(this->x_old(), x);
            // } else {
            // Create temporary for initial guess only
            auto box =
                std::make_shared<BoxConstraints<Vector_t>>(nullptr, std::make_shared<Vector_t>(obstacle_->gap()));

            auto ls = std::make_shared<LineSearchBoxProjection<Vector_t>>(box, make_ref(this->x_old()));
            alpha = ls->compute(this->x_old(), x);
            // }

            x = this->x_old() + alpha * x;

            if (alpha < 0.99) {
                utopia::out() << "initial_guess_for_solver, alpha: " << alpha << "\n";
            }

            time_derivative(x, velocity);
        }

    protected:
        const Vector_t &velocity_old() const { return velocity_old_; }

    private:
        Vector_t x_old_, velocity_old_;
        Vector_t predictor_;
        Vector_t force_old_;
        bool has_zero_density_{false};

        std::shared_ptr<IObstacle<FunctionSpace>> obstacle_;
        bool debug_{false};
        int debug_from_iteration_{0};
        bool stabilized_formulation_{true};

        std::shared_ptr<LogBarrierBase> barrier_;
        bool trivial_obstacle_{false};
        bool zero_initial_guess_{true};

        std::shared_ptr<LineSearchBoxProjection<Vector_t>> line_search_;
        bool enable_line_search_{false};

        void update_x(const Vector_t &velocity, Vector_t &x) const {
            x = predictor_;
            x += (0.5 * this->delta_time()) * (velocity - this->velocity_old());
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

                // {
                //     Scalar_t norm_B_g = norm2(barrier_temp);
                //     std::stringstream ss;
                //     ss << "norm_B_g: " << norm_B_g << "\n";
                //     x.comm().root_print(ss.str());
                // }

                g += barrier_temp;

                if (g.has_nan_or_inf()) {
                    this->~ObstacleStabilizedVelocityNewmark();
                    g.comm().barrier();
                    assert(false);
                    Utopia::Abort("barrier_gradient: NaN detected!");
                }
            }
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

                // {
                //     Scalar_t norm_B_H = norm2(temp_H);
                //     std::stringstream ss;
                //     ss << "norm_B_H: " << norm_B_H << "\n";
                //     x.comm().root_print(ss.str());
                // }

                H += temp_H;
            }
        }

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
    };

}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_STABILIZED_VELOCITY_NEWMARK_INTEGRATOR_HPP

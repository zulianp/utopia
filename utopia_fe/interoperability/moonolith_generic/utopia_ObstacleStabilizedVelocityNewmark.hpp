#ifndef UTOPIA_OBSTACLE_STABILIZED_VELOCITY_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_OBSTACLE_STABILIZED_VELOCITY_NEWMARK_INTEGRATOR_HPP

#include "utopia_ContactFactory.hpp"
#include "utopia_ContactInterface.hpp"
#include "utopia_FEModelFunction.hpp"

#include "utopia_LineSearchBoxProjection.hpp"
#include "utopia_LogBarrierFactory.hpp"
#include "utopia_ObstacleStabilizedNewmark.hpp"
#include "utopia_VelocityVariant.hpp"

#include <utility>

namespace utopia {

    // Velocity form of the contact-stabilized Newmark method (see ObstacleStabilizedNewmark), with the
    // contact constraint enforced by a barrier. The unknown is v_{n+1}, and x(v) = x_pred + dt/2 (v - v_hat_n),
    // with the projected velocity v_hat_n = (x_pred - x_n) / dt
    template <class FunctionSpace>
    class ObstacleStabilizedVelocityNewmark : public TimeDependentFunction<FunctionSpace>,
                                              public ContactDependentFunction<FunctionSpace>,
                                              public VelocityVariant<FunctionSpace>,
                                              public LSStrategy<typename Traits<FunctionSpace>::Vector> {
    public:
        using Super = utopia::TimeDependentFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Layout_t = typename Traits<FunctionSpace>::Layout;

        using Penalty = utopia::Penalty<Matrix_t, Vector_t>;

        void read(Input &in) override {
            Super::read(in);

            in.get("verbose", verbose_);
            in.get("debug", debug_);
            in.get("debug_from_iteration", debug_from_iteration_);
            in.get("stabilized_formulation", stabilized_formulation_);
            in.get("trivial_obstacle", trivial_obstacle_);
            in.get("zero_initial_guess", zero_initial_guess_);
            in.get("enable_line_search", enable_line_search_);
            in.get("enable_NaN_safe_line_search", enable_NaN_safe_line_search_);
            in.get("max_bisections", max_bisections_);
            in.get("allow_projection", allow_projection_);

            in.get("dumping", dumping_);
            in.get("consistent_initial_acceleration", consistent_initial_acceleration_);

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

            // The scaling matrix is set in update_constraints, once the obstacle has been assembled
            in.get("use_barrier_mass_scaling", use_barrier_mass_scaling_);
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

            barrier_->reset();

            auto box =
                std::make_shared<BoxConstraints<Vector_t>>(nullptr, std::make_shared<Vector_t>(obstacle_->gap()));

            barrier_->set_box_constraints(box);

            barrier_->set_selection(std::make_shared<Vector_t>(obstacle_->is_contact()));

            if (use_barrier_mass_scaling_) {
                barrier_->set_scaling_matrix(obstacle_->mass_matrix());
            }

            if (enable_line_search_) {
                // The trivial obstacle constrains x itself, otherwise the constraint is on x - x_old
                std::shared_ptr<Vector_t> offset = trivial_obstacle_ ? nullptr : make_ref(this->x_old());

                if (!line_search_) {
                    line_search_ = std::make_shared<LineSearchBoxProjection<Vector_t>>(box, offset);
                } else {
                    line_search_->set_box_constraints(box);
                    line_search_->set_offset_vector(offset);
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
            return get_alpha_aux(fun, g, velocity, correction, alpha);
        }

        bool get_alpha(LeastSquaresFunctionBase<Vector_t> &fun,
                       const Vector_t &g,
                       const Vector_t &velocity,
                       const Vector_t &correction,
                       Scalar_t &alpha) override {
            return get_alpha_aux(fun, g, velocity, correction, alpha);
        }

        template <class Fun>
        bool get_alpha_aux(Fun &fun,
                           const Vector_t &g,
                           const Vector_t &velocity,
                           const Vector_t &correction,
                           Scalar_t &alpha) {
            alpha = 1;
            bool ok = true;
            Vector_t work, x;

            if (line_search_) {
                x = velocity + correction;
                update_x(x, work);
                update_x(velocity, x);
                work -= x;
                ok = line_search_->get_alpha(fun, g, x, work, alpha);
            }

            if (enable_NaN_safe_line_search_) {
                ok = false;
                int t = 0;

                Vector_t new_grad;
                for (; t < max_bisections_; ++t) {
                    // Compute correction on velocty
                    work = velocity + alpha * correction;

                    // Convert to displacement
                    update_x(work, x);

                    // Check for failures in gradient!
                    if (this->function()->gradient(x, new_grad)) {
                        ok = true;
                        break;
                    }

                    alpha /= 2;
                }

                if (verbose_ && t > 1) {
                    utopia::out() << "reduced alpha to " << alpha << ", with " << t << " bisection(s)!\n";
                }

                if (t == max_bisections_ && !ok) {
                    this->space()->write("NaN.e", x);
                    this->~ObstacleStabilizedVelocityNewmark();
                    assert(false);
                    Utopia::Abort("Solution reached an unrecoverable state!");
                }
            }

            alpha *= dumping_;
            return ok;
        }

        void init_memory(const Layout_t & /*layout*/) override {}

        // Computes the predictor x_pred = x_n + dt v_hat_n, with the projected velocity v_hat_n
        void update_predictor() {
            assert(obstacle_);

            const Scalar_t dt = this->delta_time();

            if (stabilized_formulation_) {
                Vector_t temp;
                obstacle_->transform(velocity_old_, temp);

                // Admissible increment from x_old. The gap of the trivial obstacle refers to the reference
                // configuration, hence the part already used by x_old is removed
                Vector_t max_increment = obstacle_->gap();
                if (trivial_obstacle_) {
                    Vector_t t_x_old;
                    obstacle_->transform(x_old_, t_x_old);
                    max_increment -= t_x_old;
                }

                temp = utopia::min(temp, (1 / dt) * max_increment);

                obstacle_->inverse_transform(temp, projected_velocity_);
            } else {
                projected_velocity_ = velocity_old_;
            }

            predictor_ = x_old_;
            predictor_ += dt * projected_velocity_;
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

            if (consistent_initial_acceleration_ && !has_zero_density_ && !compute_initial_force()) {
                return false;
            }

            if (!obstacle_) {
                Utopia::Abort("ObstacleStabilizedVelocityNewmark: an obstacle or a contact is required!");
            }

            if (trivial_obstacle_) {
                // The trivial obstacle bounds the total displacement, so its gap refers to the reference configuration
                Vector_t zero(vlo, 0.);
                update_constraints(zero);
            } else {
                update_constraints(x);
            }

            update_predictor();
            return true;
        }

        // The force F(x_0) = r(x_0) - f(t_0) plays the role of -M a_0 (equation of motion at t_0)
        bool compute_initial_force() {
            if (this->has_time()) {
                this->space()->update(*this->time());
            }

            return this->function()->gradient(x_old_, force_old_);
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

            // Advances the clock (if the driver did not) so that the force is evaluated at the time of x
            bool ok = Super::update_IVP(x);

            this->function()->gradient(x, force_old_);

            // Store current solution
            x_old_ = x;
            velocity_old_ = velocity;

            if (!trivial_obstacle_) {
                update_constraints(x);
            }

            update_predictor();
            barrier_->reset();
            return ok;
        }

        bool time_derivative(const Vector_t &x, Vector_t &velocity) const override {
            const Scalar_t dt = this->delta_time();
            velocity = projected_velocity_;
            velocity += (2 / dt) * (x - predictor_);
            return true;
        }

        template <class... Args>
        ObstacleStabilizedVelocityNewmark(Args &&...args) : Super(std::forward<Args>(args)...) {}

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

        bool update(const Vector_t &velocity) override {
            if (barrier_) {
                Vector_t x, x_obs;
                update_x(velocity, x);
                constraint_coordinates(x, x_obs);
                barrier_->update(x_obs);
            }

            return true;
        }

        // Incremental potential plus barrier energy, evaluated at x(v):
        // E(x) = Psi(x) - f^T x + F_n^T x + 2/dt^2 (x - x_pred)^T M (x - x_pred) + B(x)
        bool value(const Vector_t &velocity, Scalar_t &v) const override {
            Vector_t x;
            update_x(velocity, x);

            if (!this->function()->value(x, v)) {
                return false;
            }

            if (!has_zero_density_) {
                const Scalar_t dt = this->delta_time();

                Vector_t dx = x - predictor_;
                Vector_t Mdx = (*this->mass_matrix()) * dx;
                const Scalar_t inertia = (2 / (dt * dt)) * dot(dx, Mdx);
                const Scalar_t previous_force = dot(force_old_, x);
                v += inertia + previous_force;
            }

            if (barrier_) {
                Vector_t barrier_x(layout(x), 0);
                constraint_coordinates(x, barrier_x);

                Scalar_t barrier_value = 0;
                if (!barrier_->value(barrier_x, barrier_value)) {
                    return false;
                }

                v += barrier_value;
            }

            return true;
        }

        inline Vector_t &x_old() { return x_old_; }
        inline const Vector_t &x_old() const { return x_old_; }

        const Vector_t &solution() const override { return x_old(); }
        const Vector_t &velocity() const { return velocity_old_; }

        // The displacement is written by report_solution
        bool register_output(IO<FunctionSpace> &output) override {
            set_velocity_field();
            output.register_output_field(velocity_field_);
            return true;
        }

        bool update_output(IO<FunctionSpace> &output) override {
            set_velocity_field();
            output.update_output_field(velocity_field_);
            return true;
        }

        bool report_solution(const Vector_t &) override { return Super::report_solution(solution()); }

        bool set_initial_condition(const Vector_t &x) override {
            x_old_ = x;
            return true;
        }

        inline void set_contact(const std::shared_ptr<ContactInterface<FunctionSpace>> obstacle) override {
            obstacle_ = obstacle;
        }

        // Called by Newton for every iterate. Optional (as in ObstacleVelocityNewmark): a penalty (ShiftedPenalty)
        // needs penetrating iterates to build up the contact force, and the projection would undo every step
        bool project_onto_feasibile_region(Vector_t &velocity) const final {
            bool ok = true;
            if (barrier_ && allow_projection_) {
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
            } else {
                // x(v_hat_n) is the predictor, which is admissible in the stabilized formulation
                velocity = projected_velocity_;
            }

            ensure_admissible(velocity);
            impose_velocity_constraints(velocity);
        }

        // The barrier is only defined for admissible displacements. If x(v) is not admissible (e.g., when
        // the guess moves the body beyond the gap), the guess x(v) = x_old (no motion), which always is, is used
        void ensure_admissible(Vector_t &velocity) const {
            if (!barrier_ || !barrier_->box()) return;

            Vector_t x, t_x;
            update_x(velocity, x);
            constraint_coordinates(x, t_x);

            if (barrier_->box()->count_violations(t_x) > 0) {
                time_derivative(x_old_, velocity);
            }
        }

        // Sets the Dirichlet dofs of the velocity such that x(v) = g(t_{n+1}) there,
        // i.e., v_{n+1} = v_n + 2/dt (g(t_{n+1}) - x_pred)
        void impose_velocity_constraints(Vector_t &velocity) const {
            Vector_t g = x_old_;
            this->space()->apply_constraints(g);

            Vector_t v_bc;
            time_derivative(g, v_bc);

            // 1 on the free dofs and 0 on the Dirichlet dofs
            Vector_t is_free(layout(velocity), 1.);
            this->space()->apply_zero_constraints(is_free);

            Vector_t v_free = e_mul(is_free, velocity);
            Vector_t v_bc_free = e_mul(is_free, v_bc);

            velocity = v_free;
            velocity += v_bc;
            velocity -= v_bc_free;
        }

        inline std::shared_ptr<LSStrategy<Vector_t>> line_search() override {
            if (line_search_ || enable_NaN_safe_line_search_) {
                return make_ref(*this);
            } else {
                return nullptr;
            }
        }

    protected:
        const Vector_t &velocity_old() const { return velocity_old_; }

        UTOPIA_NVCC_PRIVATE
        Vector_t x_old_, velocity_old_;
        Vector_t projected_velocity_;
        Vector_t predictor_;
        Field<FunctionSpace> velocity_field_{"velocity"};
        Vector_t force_old_;
        bool has_zero_density_{false};

        std::shared_ptr<ContactInterface<FunctionSpace>> obstacle_;
        bool debug_{false};
        int debug_from_iteration_{0};
        bool stabilized_formulation_{true};

        std::shared_ptr<Penalty> barrier_;
        bool trivial_obstacle_{false};
        bool zero_initial_guess_{true};

        std::shared_ptr<LineSearchBoxProjection<Vector_t>> line_search_;
        bool enable_line_search_{false};
        bool enable_NaN_safe_line_search_{false};
        int max_bisections_{6};
        bool verbose_{false};
        bool allow_projection_{false};

        Scalar_t dumping_{0.98};
        bool use_barrier_mass_scaling_{true};
        bool consistent_initial_acceleration_{true};

        // class Debug {
        // public:
        //     using IO_t = utopia::IO<FunctionSpace>;
        //     Debug(FunctionSpace &space) : io(std::make_shared<IO_t>(space)) {
        //         io->set_output_path("barrier_gradient.e");
        //     }

        //     void report_barrier_gradient(const Vector_t &barrier_g) { io->write(barrier_g, counter, counter++); }

        //     std::shared_ptr<IO_t> io;
        //     int counter{0};
        // };

        void update_x(const Vector_t &velocity, Vector_t &x) const {
            x = predictor_;
            x += (0.5 * this->delta_time()) * (velocity - projected_velocity_);
        }

        void set_velocity_field() {
            velocity_field_.set_space(this->space());
            velocity_field_.set_tensor_size(this->space()->n_var());
            velocity_field_.set_data(std::make_shared<Vector_t>(velocity_old_));
        }

        // Maps the displacement to the quantity bounded by the gap (x for the trivial obstacle, x - x_old otherwise)
        // in the obstacle coordinate system
        void constraint_coordinates(const Vector_t &x, Vector_t &out) const {
            if (trivial_obstacle_) {
                obstacle_->transform(x, out);
            } else {
                Vector_t delta_x = x - this->x_old();
                obstacle_->transform(delta_x, out);
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

                if (g.has_nan_or_inf()) {
                    this->space()->write("NaN_barrier_g.e", barrier_temp);
                    this->space()->write("NaN_barrier_x.e", x);

                    this->~ObstacleStabilizedVelocityNewmark();
                    g.comm().barrier();
                    assert(false);
                    Utopia::Abort("ObstacleStabilizedVelocityNewmark::barrier_gradient: NaN detected!");
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

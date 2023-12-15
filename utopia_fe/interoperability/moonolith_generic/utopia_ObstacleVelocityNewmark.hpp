#ifndef UTOPIA_OBSTACLE_VELOCITY_NEWMARK_HPP
#define UTOPIA_OBSTACLE_VELOCITY_NEWMARK_HPP

#include "utopia_LogBarrierFactory.hpp"

#include "utopia_BoxConstrainedFEFunction.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_ContactInterface.hpp"
#include "utopia_fe_Core.hpp"

#include "utopia_AnalyticObstacle_impl.hpp"
#include "utopia_ImplicitObstacle_impl.hpp"

#include "utopia_ContactFactory.hpp"

#include "utopia_IsNotSupported.hpp"
#include "utopia_LineSearchBoxProjection.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ObstacleVelocityNewmark final : public NewmarkIntegrator<FunctionSpace>,
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

        ObstacleVelocityNewmark(const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &unconstrained)
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

                // utopia::out() << "ISCONTACT\n";
                this->space()->write("isc_" + std::to_string(iter_debug) + ".e", obstacle_->is_contact());
            }
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
            bool ok = true;
            alpha = damping_;
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

                if (verbose_ && t > 1 && !g.comm().rank()) {
                    utopia::out() << "reduced alpha to " << alpha << ", with " << t << " bisection(s)!\n";
                }

                if (t == max_bisections_ && !ok) {
                    this->space()->write("NaN.e", x);
                    this->~ObstacleVelocityNewmark();
                    assert(false);
                    Utopia::Abort("Solution reached an unrecoverable state!");
                }
            }

            return ok;
        }

        void init_memory(const Layout_t & /*layout*/) override {}

        void initial_guess_for_solver(Vector_t &velocity) override {
            velocity.set(0.);

            project_onto_feasibile_region(velocity);
        }

        inline std::shared_ptr<LSStrategy<Vector_t>> line_search() override {
            if (line_search_ || enable_NaN_safe_line_search_) {
                return make_ref(*this);
            } else {
                return nullptr;
            }
        }

        void read(Input &in) override {
            Super::read(in);

            in.get("verbose", verbose_);
            in.get("debug", debug_);
            in.get("debug_from_iteration", debug_from_iteration_);
            in.get("trivial_obstacle", trivial_obstacle_);
            in.get("enable_line_search", enable_line_search_);
            in.get("enable_NaN_safe_line_search", enable_NaN_safe_line_search_);
            in.get("max_bisections", max_bisections_);
            in.get("zero_initial_guess", zero_initial_guess_);
            in.get_deprecated("dumping", "damping", damping_);
            in.get("damping", damping_);
            in.get("allow_projection", allow_projection_);
            in.get("non_smooth_projection", non_smooth_projection_);
            in.get("max_projection_iterations", max_projection_iterations_);

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

            bool use_barrier_mass_scaling = true;
            in.get("use_barrier_mass_scaling", use_barrier_mass_scaling);

            if (use_barrier_mass_scaling) {
                barrier_->set_scaling_matrix(obstacle_->mass_matrix());
            }

            if (mpi_world_rank() == 0) {
                describe(utopia::out().stream());
            }
        }

        void describe(std::ostream &os) const {
            os << "-----------------------------------------\n";
            os << "utopia::ObstacleVelocityNewmark\n";
            os << "-----------------------------------------\n";

            os << "debug:\t" << debug_ << "\n";
            os << "debug_from_iteration:\t" << debug_from_iteration_ << "\n";
            os << "trivial_obstacle:\t" << trivial_obstacle_ << "\n";
            os << "enable_line_search:\t" << enable_line_search_ << "\n";
            os << "enable_NaN_safe_line_search:\t" << enable_NaN_safe_line_search_ << "\n";
            os << "allow_projection:\t" << allow_projection_ << "\n";
            os << "max_bisections:\t" << max_bisections_ << "\n";
            os << "verbose:\t" << verbose_ << "\n";
            os << "zero_initial_guess:\t" << zero_initial_guess_ << "\n";
            os << "non_smooth_projection:\t" << non_smooth_projection_ << "\n";

            os << "damping:\t" << damping_ << "\n";

            os << "-----------------------------------------\n";
        }

        bool setup_IVP(Vector_t &x) override {
            // update_constraints(x);

            // if (!this->assemble_mass_matrix()) {
            //     return false;
            // }

            // if (non_smooth_projection_) {
            //     non_smooth_project(x);
            // }

            // assert(this->mass_matrix());
            // Scalar_t sum_mm = sum(*this->mass_matrix());
            // this->state()->has_zero_density = sum_mm == 0.0;

            // auto vlo = layout(x);

            // // x_old_.zeros(vlo);
            // this->x_old() = x;
            // this->velocity_old().zeros(vlo);
            // this->acceleration_old().zeros(vlo);

            update_constraints(x);
            return Super::setup_IVP(x);
        }

        bool time_derivative(const Vector_t &x, Vector_t &dfdt) const override {
            return Super::time_derivative(x, dfdt);
        }

        bool non_smooth_project(Vector_t &x) {
            MPRGP<Matrix_t, Vector_t> qp_solver;

            Matrix_t H;
            Vector_t buff_1(layout(x), 0), buff_2(layout(x), 0);

            // delta_x
            buff_2 = x - this->x_old();

            // Transform x into obstacle coordinate system
            obstacle_->transform(buff_2, buff_1);

            SizeType n_violations = box_->count_violations(buff_1);
            if (!n_violations) return true;

            // remove delta_x from upper_bound
            *box_->upper_bound() -= buff_1;

            qp_solver.verbose(true);
            qp_solver.set_box_constraints(*box_);

            this->function()->hessian(x, H);
            Super::integrate_hessian(x, H);
            this->space()->apply_constraints(H);

            buff_1.set(0);
            buff_2.set(0);

            Matrix_t H_c;
            obstacle_->transform(H, H_c);
            qp_solver.max_it(max_projection_iterations_);
            qp_solver.solve(H_c, buff_2, buff_1);

            Scalar_t diff_x = norm2(buff_1);

            if (!x.comm().rank()) {
                utopia::out() << "found " << n_violations << " violations, diff_x: " << diff_x << "\n";
            }

            // Transform correction into body coordinate system
            obstacle_->inverse_transform(buff_1, buff_2);
            x += buff_2;
            return true;
        }

        bool update_IVP(const Vector_t &velocity) override {
            Vector_t x = this->x_old();
            update_x(velocity, x);

            if (non_smooth_projection_) {
                non_smooth_project(x);
            }

            if (x.has_nan_or_inf()) {
                this->~ObstacleVelocityNewmark();
                x.comm().barrier();
                assert(false);
                Utopia::Abort("update_IVP: NaN detected!");
            }

            if (!trivial_obstacle_) {
                update_constraints(x);
            }

            barrier_->reset();
            return Super::update_IVP(x);
        }

        bool update_BVP() override {
            this->space()->apply_constraints(this->x_old());
            return true;
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

                if (g.has_nan_or_inf()) {
                    this->space()->write("NaN_barrier_g.e", barrier_temp);
                    this->space()->write("NaN_barrier_x.e", x);

                    this->~ObstacleVelocityNewmark();
                    g.comm().barrier();
                    assert(false);
                    Utopia::Abort("ObstacleVelocityNewmark::barrier_gradient: NaN detected!");
                }
            }
        }

        bool gradient(const Vector_t &velocity, Vector_t &g) const override {
            Vector_t x;
            update_x(velocity, x);

            if (!this->function()->gradient(x, g)) {
                return false;
            }

            Vector_t mom = velocity - this->velocity();
            mom *= 2 / this->delta_time();
            mom -= this->acceleration();

            g += (*this->mass_matrix()) * mom;

            barrier_gradient(x, g);

            if (this->must_apply_constraints_to_assembled()) {
                this->space()->apply_zero_constraints(g);
            }

            return true;
        }

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            Super::integrate_gradient(x, g);
            barrier_gradient(x, g);
        }

        bool hessian(const Vector_t &velocity, Matrix_t &H) const override {
            Vector_t x;
            update_x(velocity, x);

            if (!Super::hessian(x, H)) {
                return false;
            }

            return true;
        }

        void integrate_hessian(const Vector_t &x, Matrix_t &H) const override {
            UTOPIA_TRACE_REGION_BEGIN("ObstacleVelocityNewmark::integrate_hessian");
            Super::integrate_hessian(x, H);
            barrier_hessian(x, H);
            H *= (this->delta_time() / 2);

            UTOPIA_TRACE_REGION_END("ObstacleVelocityNewmark::integrate_hessian");
        }

        bool hessian_and_gradient(const Vector_t &velocity, Matrix_t &H, Vector_t &g) const override {
            Vector_t x;
            update_x(velocity, x);
            return Super::hessian_and_gradient(x, H, g);
        }

        ////////////////////////////////////////////////////////////////////////////////

        const Vector_t &solution() const override { return this->x_old(); }
        bool report_solution(const Vector_t &) override { return Super::report_solution(solution()); }

        bool update(const Vector_t &velocity) override {
            if (barrier_) {
                Vector_t x, x_obs;
                update_x(velocity, x);

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

        bool project_onto_feasibile_region(Vector_t &velocity) const final {
            bool ok = true;
            if (barrier_ && allow_projection_) {
                Vector_t x;
                update_x(velocity, x);

                project_x_onto_feasibile_region(x);

                Super::time_derivative(x, velocity);

                velocity.comm().root_print("project_onto_feasibile_region");
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

    private:
        std::shared_ptr<ContactInterface<FunctionSpace>> obstacle_;
        std::shared_ptr<Penalty> barrier_;
        bool debug_{false};
        int debug_from_iteration_{0};
        bool trivial_obstacle_{false};
        bool enable_line_search_{false};
        bool enable_NaN_safe_line_search_{false};
        bool allow_projection_{false};
        int max_bisections_{6};
        bool verbose_{false};
        bool zero_initial_guess_{true};
        bool non_smooth_projection_{false};
        int max_projection_iterations_{10000};

        Scalar_t damping_{1};

        std::shared_ptr<LineSearchBoxProjection<Vector_t>> line_search_;
        std::shared_ptr<BoxConstraints<Vector_t>> box_;

        void update_x(const Vector_t &velocity, Vector_t &x) const {
            x = this->x_old();
            x += (0.5 * this->delta_time()) * (this->velocity_old() + velocity);
        }

        bool update_constraints(const Vector_t &x) {
            x.comm().root_print("ObstacleVelocityNewmark::update_constraints\n");

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
            barrier_->set_scaling_matrix(obstacle_->mass_matrix());

            if (enable_line_search_) {
                if (!line_search_) {
                    line_search_ = std::make_shared<LineSearchBoxProjection<Vector_t>>(box, make_ref(this->x_old()));
                } else {
                    line_search_->set_box_constraints(box);

                    if (!trivial_obstacle_) {
                        line_search_->set_offset_vector(make_ref(this->x_old()));
                    }
                }

                auto trafo = obstacle_->orthogonal_transformation();
                assert(trafo);
                if (!trafo) {
                    Utopia::Abort(
                        "ObstacleVelocityNewmark:update_constraints: orthogonal_transformation is mandatory for "
                        "line_search!");
                }

                line_search_->set_transform(trafo);
                line_search_->set_dumping(damping_);
            }

            box_ = box;

            return ok;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_VELOCITY_NEWMARK_HPP

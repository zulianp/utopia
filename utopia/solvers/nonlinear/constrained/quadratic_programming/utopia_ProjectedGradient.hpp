#ifndef UTOPIA_PROJECTED_GRADIENT_HPP
#define UTOPIA_PROJECTED_GRADIENT_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Layout.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_Operations.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_Recorder.hpp"
#include "utopia_Smoother.hpp"

#include <cassert>
#include <cmath>

namespace utopia {
    // slow and innefficient implementation just for testing
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class ProjectedGradient final : public OperatorBasedQPSolver<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using Super = utopia::OperatorBasedQPSolver<Matrix, Vector>;

        using Super::solve;
        using Super::update;

        ProjectedGradient() = default;

        ProjectedGradient(const ProjectedGradient &) = default;

        inline ProjectedGradient *clone() const override {
            auto ptr = new ProjectedGradient(*this);
            ptr->set_box_constraints(this->get_box_constraints());

            return ptr;
        }

        void read(Input &in) override { OperatorBasedQPSolver<Matrix, Vector>::read(in); }

        void print_usage(std::ostream &os) const override { OperatorBasedQPSolver<Matrix, Vector>::print_usage(os); }

        bool apply(const Vector &b, Vector &x) override {
            return solve(operator_cast<Vector>(*this->get_operator()), b, x);
        }

        bool solve_unconstrained(const Operator<Vector> &A, const Vector &b, Vector &x) {
            // UTOPIA_RECORD_SCOPE_BEGIN("apply");

            if (this->verbose()) this->init_solver("utopia ProjectedGradient", {" it. ", "|| u - u_old ||"});

            init_memory(layout(b));

            x_old = x;
            A.apply(x, u);
            p = b - u;
            Scalar alpha = 1.;

            bool converged = false;
            const SizeType check_s_norm_each = 20;

            int iteration = 0;
            while (!converged) {
                // perform step
                x = x + alpha * p;

                A.apply(x, u);
                p = b - u;

                if (iteration % check_s_norm_each == 0) {
                    const Scalar diff = norm2(x_old - x);

                    if (this->verbose()) {
                        PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
                    }

                    converged = this->check_convergence(iteration, 1, 1, diff);
                }

                ++iteration;

                if (converged) break;

                x_old = x;
                A.apply(p, Ap);
                alpha = dot(p, p) / dot(p, Ap);

                if (std::isinf(alpha) || alpha == 0. || std::isnan(alpha)) {
                    const Scalar diff = norm2(x_old - x);

                    // UTOPIA_RECORD_VALUE("x_old - x", Vector(x_old - x));

                    if (this->verbose()) {
                        PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
                    }

                    converged = this->check_convergence(iteration, 1, 1, diff);
                    break;
                }
            }
            // UTOPIA_RECORD_SCOPE_END("apply");
            return converged;
        }

        bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override {
            // UTOPIA_RECORD_SCOPE_BEGIN("apply");

            if (!this->has_bound()) {
                return solve_unconstrained(A, b, x);
            }

            if (this->verbose()) this->init_solver("utopia ProjectedGradient", {" it. ", "|| u - u_old ||"});

            init_memory(layout(b));

            // ideally, we have two separate implementations, or cases
            this->fill_empty_bounds(layout(x));

            const auto &upbo = this->get_upper_bound();
            const auto &lobo = this->get_lower_bound();

            x_old = x;
            A.apply(x, u);
            u = b - u;
            // u = b - A * x;
            p = u;
            Scalar alpha = 1.;

            // UTOPIA_RECORD_VALUE("u = b - A * x", u);
            // UTOPIA_RECORD_VALUE("p = u", p);

            bool converged = false;
            const SizeType check_s_norm_each = 20;

            int iteration = 0;
            while (!converged) {
                // perform step
                x_half = x + alpha * p;

                // UTOPIA_RECORD_VALUE("x_half = x + alpha * p", x_half);

                x = utopia::max(lobo, x_half);

                // UTOPIA_RECORD_VALUE("x = utopia::max(lobo, x_half)", x);

                x = utopia::min(upbo, x);

                // UTOPIA_RECORD_VALUE("x = utopia::min(upbo, x)", x);

                A.apply(x, u);
                u = b - u;

                // u = b - A * x;
                {
                    Read<Vector> r_u(u), r_x(x), r_upbo(upbo), r_lobo(lobo);
                    Write<Vector> w_(p, utopia::LOCAL);

                    auto r = range(x);

                    for (auto i = r.begin(); i != r.end(); ++i) {
                        const auto x_i = x.get(i);

                        if (approxeq(x_i, upbo.get(i)) || approxeq(x_i, lobo.get(i))) {
                            p.set(i, 0);
                        } else {
                            p.set(i, u.get(i));
                        }
                    }
                }

                // UTOPIA_RECORD_VALUE("p <- min_max", p);

                if (iteration % check_s_norm_each == 0) {
                    const Scalar diff = norm2(x_old - x);

                    // UTOPIA_RECORD_VALUE("x_old - x", Vector(x_old - x));

                    if (this->verbose()) {
                        PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
                    }

                    converged = this->check_convergence(iteration, 1, 1, diff);
                }

                ++iteration;

                if (converged) break;

                x_old = x;
                A.apply(p, Ap);
                alpha = dot(u, p) / dot(p, Ap);

                if (std::isinf(alpha) || alpha == 0. || std::isnan(alpha)) {
                    const Scalar diff = norm2(x_old - x);

                    // UTOPIA_RECORD_VALUE("x_old - x", Vector(x_old - x));

                    if (this->verbose()) {
                        PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
                    }

                    converged = this->check_convergence(iteration, 1, 1, diff);
                    break;
                }
            }
            // UTOPIA_RECORD_SCOPE_END("apply");
            return converged;
        }

        void init_memory(const Layout &layout) override {
            OperatorBasedQPSolver<Matrix, Vector>::init_memory(layout);

            p.zeros(layout);
            Ap.zeros(layout);
            x_old.zeros(layout);
            x_half.zeros(layout);
        }

        void update(const Operator<Vector> &) override {}

    private:
        // buffers
        Vector x_old, x_half, p, u, Ap;
    };
}  // namespace utopia

#endif  // UTOPIA_PROJECTED_GRADIENT_HPP

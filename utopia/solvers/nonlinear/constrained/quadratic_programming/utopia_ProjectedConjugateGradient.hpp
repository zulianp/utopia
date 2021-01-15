#ifndef UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP
#define UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP

#include "utopia_ForwardDeclarations.hpp"

#include "utopia_Algorithms.hpp"

#include <cassert>
#include <cmath>

namespace utopia {
    // slow and innefficient implementation just for testing
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class ProjectedConjugateGradient : public QPSolver<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        ProjectedConjugateGradient() = default;

        ProjectedConjugateGradient(const ProjectedConjugateGradient &other)
            : VariableBoundSolverInterface<Vector>(other),
              PreconditionedSolverInterface<Vector>(other),
              QPSolver<Matrix, Vector>(other) {}

        inline ProjectedConjugateGradient *clone() const override {
            auto ptr = new ProjectedConjugateGradient(*this);
            ptr->set_box_constraints(this->get_box_constraints());
            return ptr;
        }

        bool apply(const Vector &b, Vector &x) override {
            if (this->verbose()) {
                this->init_solver("ProjectedConjugateGradient", {" it. ", "|| u - u_old ||"});
            }

            const Matrix &A = *this->get_operator();

            // ideally, we have two separate implementations, or cases
            this->fill_empty_bounds(layout(x));

            const auto &ub = this->get_upper_bound();
            const auto &lb = this->get_lower_bound();

            x_old = x;
            uk = b - A * x;

            bool converged = false;
            const SizeType check_s_norm_each = 10;
            pk = -uk;

            int iteration = 0;
            while (!converged) {
                // START step

                Apk = A * pk;
                Scalar alpha = dot(uk, pk) / dot(pk, Apk);
                assert(alpha != 0.);
                if (alpha == 0. || std::isinf(alpha) || std::isnan(alpha)) break;

                x_half = x_old + alpha * pk;

                x = utopia::min(x_half, ub);
                x = utopia::max(x, lb);

                uk = b - A * x;

                {
                    auto uk_view = const_local_view_device(uk);
                    auto ub_view = const_local_view_device(ub);
                    auto lb_view = const_local_view_device(lb);
                    auto pk_view = const_local_view_device(pk);

                    auto x_view = local_view_device(zk);
                    auto wk_view = local_view_device(wk);
                    auto zk_view = local_view_device(zk);

                    parallel_for(
                        local_range_device(x), UTOPIA_LAMBDA(const SizeType &i) {
                            const auto elem = x_view.get(i);

                            Scalar val = 0.;
                            if (device::approxeq(elem, ub_view.get(i), device::epsilon<Scalar>()) ||
                                device::approxeq(elem, lb_view.get(i), device::epsilon<Scalar>())) {
                                val = device::max(uk_view.get(i), Scalar(0));
                            } else {
                                val = uk_view.get(i);
                            }

                            if (val == 0) {
                                zk_view.set(i, device::max(pk_view.get(i), Scalar(0)));
                            } else {
                                zk_view.set(i, pk_view.get(i));
                            }

                            wk_view.set(i, val);
                        });
                }

                Apk = A * pk;
                const Scalar beta = dot(wk, A * pk) / dot(pk, Apk);
                pk = wk + beta * zk;

                // END step

                if (iteration % check_s_norm_each == 0 || std::isinf(beta) || std::isnan(beta)) {
                    const Scalar diff = norm2(x_old - x);

                    if (this->verbose()) {
                        PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
                    }

                    converged = this->check_convergence(iteration, 1, 1, diff);
                }

                ++iteration;

                if (converged) break;

                x_old = x;
            }

            return converged;
        }

        void init_memory(const Layout &layout) override {
            QPSolver<Matrix, Vector>::init_memory(layout);
            r.zeros(layout);
            uk.zeros(layout);
            wk.zeros(layout);
            zk.zeros(layout);
            pk.zeros(layout);
            Apk.zeros(layout);
        }

        void update(const std::shared_ptr<const Matrix> &op) override {
            QPSolver<Matrix, Vector>::update(op);
            init_memory(row_layout(*op));
        }

    private:
        // buffers
        Vector x_old, x_half, r, uk, wk, zk, pk, Apk;
    };
}  // namespace utopia

#endif  // UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP

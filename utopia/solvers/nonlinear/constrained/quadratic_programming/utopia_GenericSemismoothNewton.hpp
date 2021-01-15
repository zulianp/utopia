#ifndef UTOPIA_GENERIC_SEMISMOOTH_NEWTON_HPP
#define UTOPIA_GENERIC_SEMISMOOTH_NEWTON_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"

#include "utopia_Core.hpp"
#include "utopia_Wrapper.hpp"

#include <memory>
#include <vector>

namespace utopia {

    template <class Matrix, class Vector, class Constraint, int Backend = Traits<Vector>::Backend>
    class GenericSemismoothNewton : public IterativeSolver<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> Solver;

        GenericSemismoothNewton(const Constraint &constraint, const std::shared_ptr<Solver> &linear_solver)
            : constraint_(constraint), linear_solver_(linear_solver) {}

        GenericSemismoothNewton *clone() const override {
            return new GenericSemismoothNewton(constraint_, std::shared_ptr<Solver>(linear_solver_->clone()));
        }

        void read(Input &in) override {
            IterativeSolver<Matrix, Vector>::read(in);

            if (linear_solver_) {
                in.get("linear-solver", *linear_solver_);
            }
        }

        void print_usage(std::ostream &os) const override {
            IterativeSolver<Matrix, Vector>::print_usage(os);
            this->print_param_usage(os, "linear-solver", "LinearSolver", "Input parameters for linear solver.", "-");
        }

        bool solve(const Matrix &A, const Vector &b, Vector &x_new) override {
            // const SizeType local_N = local_size(x_new).get(0);

            auto vec_layout = layout(b);
            auto mat_layout = square_matrix_layout(vec_layout);

            SizeType iterations = 0;
            bool converged = false;

            Vector active(vec_layout, 0.0);
            Vector g(vec_layout, 0.0);
            Vector prev_active(vec_layout, 0.0);
            Vector x_old = x_new;

            // active/inactive constraints
            Matrix A_c;
            Matrix I_c;

            if (is_sparse<Matrix>::value) {
                A_c.sparse(mat_layout, 1, 0.0);
                I_c.sparse(mat_layout, 1, 0.0);
            } else {
                A_c.dense(mat_layout, 0.0);
                I_c.dense(mat_layout, 0.0);
            }

            Scalar f_norm = 9e9;

            if (this->verbose()) this->init_solver("SEMISMOOTH NEWTON METHOD", {" it. ", "|| g ||"});

            while (!converged) {
                // reminder: complementarity-condition
                constraint_(A, b, x_new, active, g);

                {
                    Write<Matrix> w_A_c(A_c);
                    Write<Matrix> w_I_c(I_c);
                    Read<Vector> r_v(active);

                    const Range rr = row_range(A_c);

                    for (SizeType i = rr.begin(); i != rr.end(); i++) {
                        if (active.get(i) > 1e-8) {
                            A_c.set(i, i, 1.0);
                            I_c.set(i, i, 0.0);
                        } else {
                            I_c.set(i, i, 1.0);
                            A_c.set(i, i, 0.0);
                        }
                    }
                }

                if (iterations > 0) {
                    const SizeType n_changed = Scalar(norm1(prev_active - active));

                    if (n_changed == 0) {
                        // active set doesn't change anymore => converged
                        // fix this to be done in other way
                        converged = this->check_convergence(iterations, 1e-15, 1, 1);
                        return true;
                    }
                }

                prev_active = active;

                Matrix H = A_c + I_c * A;
                Vector sub_g = (I_c * b + A_c * g);

                assert(!has_nan_or_inf(H));
                assert(!has_nan_or_inf(g));

                if (!linear_solver_->solve(H, sub_g, x_new)) {
                    std::cerr << "[Error] linear solver did not manage to solve the linear system" << std::endl;
                    break;
                }

                if (has_nan_or_inf(x_new)) {
                    assert(!has_nan_or_inf(x_new));
                    std::cerr << "[Error] nan/inf entries in the solution vector" << std::endl;
                    converged = false;
                    break;
                }

                f_norm = norm2(x_new - x_old);

                // print iteration status on every iteration
                if (this->verbose()) {
                    PrintInfo::print_iter_status(iterations, {f_norm});
                }

                converged = this->check_convergence(iterations, f_norm, 1, 1);

                x_old = x_new;
                iterations++;
            }

            return converged;
        }

        Constraint constraint_;
        std::shared_ptr<Solver> linear_solver_;
    };
}  // namespace utopia

#endif  // UTOPIA_GENERIC_SEMISMOOTH_NEWTON_HPP

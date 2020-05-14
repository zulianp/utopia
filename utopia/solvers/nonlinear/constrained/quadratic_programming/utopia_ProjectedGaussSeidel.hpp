#ifndef UTOPIA_PROJECTED_GAUSS_SEIDEL_HPP
#define UTOPIA_PROJECTED_GAUSS_SEIDEL_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Layout.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_RowView.hpp"
#include "utopia_Smoother.hpp"

#include <cmath>

namespace utopia {
    // slow and innefficient implementation just for testing
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class ProjectedGaussSeidel : public QPSolver<Matrix, Vector> {
       public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        ProjectedGaussSeidel() : n_local_sweeps_(3) {}

        ProjectedGaussSeidel(const ProjectedGaussSeidel &) = default;

        inline ProjectedGaussSeidel *clone() const override {
            auto ptr = new ProjectedGaussSeidel(*this);
            ptr->set_box_constraints(this->get_box_constraints());
            return ptr;
        }

        void read(Input &in) override {
            QPSolver<Matrix, Vector>::read(in);

            in.get("use_line_search", use_line_search_);
            in.get("use_symmetric_sweep", use_symmetric_sweep_);
            in.get("n_local_sweeps", n_local_sweeps_);
            in.get("l1", l1_);
        }

        void print_usage(std::ostream &os) const override {
            QPSolver<Matrix, Vector>::print_usage(os);

            this->print_param_usage(os, "use_line_search", "bool", "Determines if line-search should be used.", "true");
            this->print_param_usage(
                os, "use_symmetric_sweep", "bool", "Determines if symmetric local should be used.", "true");
            this->print_param_usage(os, "n_local_sweeps", "int", "Number of local sweeps.", "3");
        }

        bool smooth(const Vector &b, Vector &x) override {
            const Matrix &A = *this->get_operator();

            // init(A);
            SizeType it = 0;
            SizeType n_sweeps = this->sweeps();
            if (this->has_bound()) {
                while (it++ < n_sweeps) {
                    step(A, b, x);
                    std::cout
                        << "--------------------------------------------- sweep --------------------------------- \n";
                }
            } else {
                while (unconstrained_step(A, b, x) && it++ < n_sweeps) {
                }
            }
            return it == SizeType(this->sweeps() - 1);
        }

        bool apply(const Vector &b, Vector &x) override {
            if (this->verbose()) {
                if (l1_) {
                    this->init_solver("utopia ProjectedL1GaussSeidel", {" it. ", "|| u - u_old ||"});
                } else {
                    this->init_solver("utopia ProjectedGaussSeidel", {" it. ", "|| u - u_old ||"});
                }
            }

            const Matrix &A = *this->get_operator();

            x_old = x;
            bool converged = false;
            const SizeType check_s_norm_each = 1;

            int iteration = 0;
            while (!converged) {
                if (this->has_bound()) {
                    step(A, b, x);
                } else {
                    unconstrained_step(A, b, x);
                }

                if (iteration % check_s_norm_each == 0) {
                    c = x - x_old;
                    const Scalar diff = norm2(c);

                    if (this->verbose()) {
                        PrintInfo::print_iter_status(iteration, {diff});
                    }

                    converged = this->check_convergence(iteration, 1, 1, diff);
                }

                ++iteration;

                if (converged) break;

                x_old = x;
            }
            return converged;
        }

        void non_linear_jacobi_step(const Matrix &A, const Vector &b, Vector &x) {
            r = b - A * x;
            x = min(x + e_mul(d_inv, r), this->get_upper_bound());
        }

        bool unconstrained_step(const Matrix &A, const Vector &b, Vector &x) {
            r = b - A * x;
            c *= 0.;

            Range rr = row_range(A);
            {
                ReadAndWrite<Vector> rw_c(c);
                Read<Vector> r_r(r), r_d_inv(d_inv);
                Read<Matrix> r_A(A);

                for (SizeType il = 0; il < this->n_local_sweeps(); ++il) {
                    for (auto i = rr.begin(); i != rr.end(); ++i) {
                        RowView<const Matrix> row_view(A, i);
                        decltype(i) n_values = row_view.n_values();

                        auto s = r.get(i);

                        for (auto index = 0; index < n_values; ++index) {
                            const decltype(i) j = row_view.col(index);
                            const auto a_ij = row_view.get(index);

                            if (rr.inside(j) && i != j) {
                                s -= a_ij * c.get(j);
                            }
                        }

                        // update correction
                        c.set(i, d_inv.get(i) * s);
                    }

                    if (use_symmetric_sweep_) {
                        for (auto i = rr.end() - 1; i >= rr.begin(); --i) {
                            RowView<const Matrix> row_view(A, i);
                            decltype(i) n_values = row_view.n_values();

                            auto s = r.get(i);

                            for (auto index = 0; index < n_values; ++index) {
                                const decltype(i) j = row_view.col(index);
                                const auto a_ij = row_view.get(index);

                                if (rr.inside(j) && i != j) {
                                    s -= a_ij * c.get(j);
                                }
                            }

                            // update correction
                            c.set(i, d_inv.get(i) * s);
                        }
                    }
                }
            }

            Scalar alpha = 1.;

            if (use_line_search_) {
                Ac = A * c;
                alpha = dot(c, r) / dot(Ac, c);

                if (std::isinf(alpha)) {
                    return true;
                }

                if (std::isnan(alpha)) {
                    return false;
                }

                if (alpha <= 0) {
                    std::cerr << "[Warning] negative alpha" << std::endl;
                    alpha = 1.;
                    c = r;
                }
            }

            x += alpha * c;
            return true;
        }

        virtual bool step(const Matrix &A, const Vector &b, Vector &x) {
            r = b - A * x;

            // localize gap function for correction
            this->fill_empty_bounds(layout(x));
            ub_loc = this->get_upper_bound() - x;
            lb_loc = this->get_lower_bound() - x;

            c *= 0.;

            Range rr = row_range(A);
            {
                ReadAndWrite<Vector> rw_c(c);
                Read<Vector> r_r(r), r_d_inv(d_inv), r_g(ub_loc), rl(lb_loc);
                Read<Matrix> r_A(A);

                for (SizeType il = 0; il < this->n_local_sweeps(); ++il) {
                    for (auto i = rr.begin(); i != rr.end(); ++i) {
                        RowView<const Matrix> row_view(A, i);
                        decltype(i) n_values = row_view.n_values();

                        auto s = r.get(i);

                        for (auto index = 0; index < n_values; ++index) {
                            const decltype(i) j = row_view.col(index);
                            const auto a_ij = row_view.get(index);

                            if (rr.inside(j) && i != j) {
                                s -= a_ij * c.get(j);
                            }
                        }

                        // update correction
                        c.set(i, std::max(std::min(d_inv.get(i) * s, ub_loc.get(i)), lb_loc.get(i)));
                    }

                    if (use_symmetric_sweep_) {
                        for (auto i = rr.end() - 1; i >= rr.begin(); --i) {
                            RowView<const Matrix> row_view(A, i);
                            decltype(i) n_values = row_view.n_values();

                            auto s = r.get(i);

                            for (auto index = 0; index < n_values; ++index) {
                                const decltype(i) j = row_view.col(index);
                                const auto a_ij = row_view.get(index);

                                if (rr.inside(j) && i != j) {
                                    s -= a_ij * c.get(j);
                                }
                            }

                            // update correction
                            // c.set(i, std::min( d_inv.get(i) * s, ub_loc.get(i)) );
                            c.set(i, std::max(std::min(d_inv.get(i) * s, ub_loc.get(i)), lb_loc.get(i)));
                        }
                    }
                }
            }

            if (use_line_search_) {
                UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel21");
                inactive_set_ *= 0.;
                UTOPIA_NO_ALLOC_END();

                UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel22");
                {
                    Read<Vector> r_c(c), r_g(ub_loc), r_l(lb_loc);
                    Write<Vector> w_a(inactive_set_);

                    for (auto i = rr.begin(); i != rr.end(); ++i) {
                        if ((c.get(i) < ub_loc.get(i)) && (c.get(i) > lb_loc.get(i))) {
                            inactive_set_.set(i, 1.);
                        }
                    }
                }
                UTOPIA_NO_ALLOC_END();

                UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel23");
                is_c_ = e_mul(c, inactive_set_);

                Ac = A * is_c_;
                Scalar alpha = dot(is_c_, r) / dot(Ac, is_c_);
                UTOPIA_NO_ALLOC_END();

                if (std::isinf(alpha)) {
                    return true;
                }

                if (std::isnan(alpha)) {
                    return false;
                }

                assert(alpha > 0);

                if (alpha <= 0) {
                    std::cerr << "[Warning] negative alpha" << std::endl;
                    alpha = 1.;
                    UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel24");
                    descent_dir = utopia::max(utopia::min(r, ub_loc), ub_loc);
                    UTOPIA_NO_ALLOC_END();
                } else if (alpha <= 1.) {
                    UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel25");
                    descent_dir = alpha * c;
                    UTOPIA_NO_ALLOC_END();
                } else {
                    UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel26");
                    c *= alpha;
                    descent_dir = utopia::max(utopia::min(c, ub_loc), ub_loc);
                    UTOPIA_NO_ALLOC_END();
                }
                UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel27");
                x += descent_dir;
                UTOPIA_NO_ALLOC_END();
            } else {
                UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel3");
                x += c;
                UTOPIA_NO_ALLOC_END();
            }

            return true;
        }

        void init_memory(const Layout &layout) override {
            d.zeros(layout);
            c.zeros(layout);

            if (use_line_search_) {
                inactive_set_.zeros(layout);
                Ac.zeros(layout);
                is_c_.zeros(layout);
                descent_dir.zeros(layout);
            }
        }

        void init(const Matrix &A) {
            auto A_layout = row_layout(A);

            if (empty(d) || !A_layout.same(layout(d))) {
                init_memory(A_layout);
            } else {
                c.set(0.0);

                if (use_line_search_) {
                    inactive_set_.set(0);
                }
            }

            d = diag(A);

            if (l1_) {
                Write<Vector> w(d);
                each_read(
                    A, [this](const SizeType &i, const SizeType &, const Scalar &value) { d.add(i, std::abs(value)); });
            }

            d_inv = 1. / d;
        }

        void update(const std::shared_ptr<const Matrix> &op) override {
            IterativeSolver<Matrix, Vector>::update(op);
            init(*op);
        }

        void use_line_search(const bool val) { use_line_search_ = val; }

        inline SizeType n_local_sweeps() const { return n_local_sweeps_; }

        inline void n_local_sweeps(const SizeType n_local_sweeps) { n_local_sweeps_ = n_local_sweeps; }

        inline void use_symmetric_sweep(const bool use_symmetric_sweep) { use_symmetric_sweep_ = use_symmetric_sweep; }

        inline void l1(const bool val) { l1_ = val; }

       private:
           bool use_line_search_{false};
           bool use_symmetric_sweep_{true};
           bool l1_{false};
           SizeType n_local_sweeps_;

           Vector r, d, ub_loc, lb_loc, c, d_inv, x_old, descent_dir, Ac;
           Vector inactive_set_;
           Vector is_c_;
    };
}  // namespace utopia

#endif  // UTOPIA_PROJECTED_GAUSS_SEIDEL_HPP

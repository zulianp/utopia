#ifndef UTOPIA_QUAD_CUB_BACKTRACKING_HPP
#define UTOPIA_QUAD_CUB_BACKTRACKING_HPP

#include "utopia_Function.hpp"
#include "utopia_FunctionNormalEq.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_PrintInfo.hpp"

namespace utopia {

    /**
     * @brief

        This is the algorithm described in the book
        Dennis and Schnabel "Numerical Methods for Nonlinear Equations
        and Unconstrained Optimization", Prentice-Hall (1983),
        reprinte by SIAM (1996), Section 6.3.2.
        @todo check params naming properly...
     */
    template <class Vector, int Backend = Traits<Vector>::Backend>
    class Backtracking final : public LSStrategy<Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

    public:
        Backtracking() : LSStrategy<Vector>() { this->c2(1e-8); }

        void read(Input &in) override { LSStrategy<Vector>::read(in); }

        void print_usage(std::ostream &os) const override { LSStrategy<Vector>::print_usage(os); }

        /**
         * @brief      Get the alpha_k on given iterate. We are using quadratic and
         * qubic interpolation as part of backtracking. For checking decrease
         * conditions, we are using Wolfe conditions.
         *
         * @param      fun      The fun with eval.
         * @param[in]  g        The gradient.
         * @param[in]  x        The current iterate.
         * @param[in]  d        The descent direction/ usually Newton step.
         * @param      alpha_k  The new step size
         *
         * @return
         */

        bool get_alpha(LeastSquaresFunctionBase<Vector> &fun,
                       const Vector &g,
                       const Vector &x,
                       const Vector &d,
                       Scalar &alpha) override {
            return get_alpha_aux_home_made(fun, g, x, d, alpha);
        }

        bool get_alpha(FunctionBase<Vector> &fun,
                       const Vector &g,
                       const Vector &x,
                       const Vector &d,
                       Scalar &alpha) override {
            return get_alpha_aux_home_made(fun, g, x, d, alpha);
        }


        bool get_alpha(FunctionBase<Vector> &fun, const Vector &g, const Vector &x,
                       const Vector &d, const Scalar &f_old, Scalar &f_new,
                       Vector &g_new, Scalar &alpha) {
          Scalar alpha_c, alpha_p, dg = dot(d, g);
          Scalar f, f0, fc, fp, t1, t2, t3, a, b, disc;
          alpha = 1.0;

          if (dg >= 0.0) {
            if (mpi_world_rank() == 0) {
              std::cerr << "utopia::LS::backtracking:: d is not descent direction \n";
            }

            alpha = 0.0;
            return false;
          }

          UTOPIA_NO_ALLOC_BEGIN("Backtracking 1");
          f = f_old;
          f0 = f;
          fc = f;
          alpha_c = alpha;
          UTOPIA_NO_ALLOC_END();

          Scalar it = 0;

          if (this->verbose())
            PrintInfo::print_init("BACKTRACKING_LS_INNER_ITERATIONS",
                                  {" it. ", "|| alpha ||"});

          while (alpha > this->c2() && it < this->max_it()) {
            UTOPIA_NO_ALLOC_BEGIN("Backtracking 2");
            x_k = x + alpha * d;
            fun.value_and_gradient(x_k, f, g_new);
            UTOPIA_NO_ALLOC_END();

            // check decrease condition (wolfe condition)
            if (f < f0 + this->c1() * alpha * dg) {
              f_new = f;
              return true;
            }

            UTOPIA_NO_ALLOC_BEGIN("Backtracking 3");
            alpha_p = alpha_c;
            alpha_c = alpha;
            fp = fc;
            fc = f;
            UTOPIA_NO_ALLOC_END();

            //  compute next step size alpha
            if (it == 0) {
              UTOPIA_NO_ALLOC_BEGIN("Backtracking 4");
              alpha = -dg / (2 * (fc - f0 - dg));
              it++;
              UTOPIA_NO_ALLOC_END();
            } else {
              UTOPIA_NO_ALLOC_BEGIN("Backtracking 5");
              // all subsequent backtracks: cubic fit
              t1 = fc - f0 - alpha_c * dg;
              t2 = fp - f0 - alpha_p * dg;
              t3 = 1 / (alpha_c - alpha_p);

              a = t3 * (t1 / std::pow(alpha_c, 2) - t2 / std::pow(alpha_p, 2));
              b = t3 * (t2 * alpha_c / std::pow(alpha_p, 2) -
                        t1 * alpha_p / std::pow(alpha_c, 2));
              disc = std::pow(b, 2) - 3 * a * dg;
              UTOPIA_NO_ALLOC_END();

              if (a != 0) {
                // cubic has unique minimum
                alpha = (-b + std::sqrt(disc)) / (3 * a);
              } else {
                // cubic is a quadratic
                alpha = -dg / (2 * b);
              }
            }

            UTOPIA_NO_ALLOC_BEGIN("Backtracking 6");
            //  saveguard the step size
            if (alpha > 0.5 * alpha_c) {
              alpha = 0.5 * alpha_c;
            }

            if (alpha < alpha_c / 10) {
              alpha = alpha_c / 10;
            }
            UTOPIA_NO_ALLOC_END();

            it++;
            if (this->verbose())
              PrintInfo::print_iter_status({it, alpha});
          }

          f_new = f;
          return true;
        }


    private:
        template <class FunctionT>
        bool get_alpha_aux_home_made(FunctionT &fun, const Vector &g, const Vector &x, const Vector &d, Scalar &alpha) {
            Scalar alpha_c, alpha_p, dg = dot(d, g);
            Scalar f, f0, fc, fp, t1, t2, t3, a, b, disc;
            alpha = 1.0;

            if (dg >= 0.0) {
                if (mpi_world_rank() == 0) {
                    std::cerr << "utopia::LS::backtracking:: d is not descent direction \n";
                }

                alpha = 0.0;
                return false;
            }

            UTOPIA_NO_ALLOC_BEGIN("Backtracking 1");
            fun.value(x, f);
            f0 = f;
            fc = f;
            alpha_c = alpha;
            UTOPIA_NO_ALLOC_END();

            Scalar it = 0;

            if (this->verbose()) PrintInfo::print_init("BACKTRACKING_LS_INNER_ITERATIONS", {" it. ", "|| alpha ||"});

            while (alpha > this->c2() && it < this->max_it()) {
                UTOPIA_NO_ALLOC_BEGIN("Backtracking 2");
                x_k = x + alpha * d;
                fun.value(x_k, f);
                UTOPIA_NO_ALLOC_END();

                // check decrease condition (wolfe condition)
                if (f < f0 + this->c1() * alpha * dg) {
                    return true;
                }

                UTOPIA_NO_ALLOC_BEGIN("Backtracking 3");
                alpha_p = alpha_c;
                alpha_c = alpha;
                fp = fc;
                fc = f;
                UTOPIA_NO_ALLOC_END();

                //  compute next step size alpha
                if (it == 0) {
                    UTOPIA_NO_ALLOC_BEGIN("Backtracking 4");
                    alpha = -dg / (2 * (fc - f0 - dg));
                    it++;
                    UTOPIA_NO_ALLOC_END();
                } else {
                    UTOPIA_NO_ALLOC_BEGIN("Backtracking 5");
                    // all subsequent backtracks: cubic fit
                    t1 = fc - f0 - alpha_c * dg;
                    t2 = fp - f0 - alpha_p * dg;
                    t3 = 1 / (alpha_c - alpha_p);

                    a = t3 * (t1 / std::pow(alpha_c, 2) - t2 / std::pow(alpha_p, 2));
                    b = t3 * (t2 * alpha_c / std::pow(alpha_p, 2) - t1 * alpha_p / std::pow(alpha_c, 2));
                    disc = std::pow(b, 2) - 3 * a * dg;
                    UTOPIA_NO_ALLOC_END();

                    if (a != 0) {
                        // cubic has unique minimum
                        alpha = (-b + std::sqrt(disc)) / (3 * a);
                    } else {
                        // cubic is a quadratic
                        alpha = -dg / (2 * b);
                    }
                }

                UTOPIA_NO_ALLOC_BEGIN("Backtracking 6");
                //  saveguard the step size
                if (alpha > 0.5 * alpha_c) {
                    alpha = 0.5 * alpha_c;
                }

                if (alpha < alpha_c / 10) {
                    alpha = alpha_c / 10;
                }
                UTOPIA_NO_ALLOC_END();

                it++;
                if (this->verbose()) PrintInfo::print_iter_status({it, alpha});
            }

            return true;
        }

    public:
        void init_memory(const Layout &layout) override {
            if (empty(x_k)) {
                x_k.zeros(layout);
            } else if (!x_k.comm().conjunction(layout.local_size() == x_k.local_size())) {
                x_k.zeros(layout);
            }
        }

    private:
        Vector x_k;
    };
}  // namespace utopia

#endif  // UTOPIA_QUAD_CUB_BACKTRACKING_HPP

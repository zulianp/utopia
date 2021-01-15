#ifndef UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP
#define UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP

#include "utopia_CauchyPoint.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_TRSubproblem.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class Dogleg final : public TRSubproblem<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

    public:
        Dogleg(
            const std::shared_ptr<LinearSolver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >())
            : TRSubproblem<Matrix, Vector>(), ls_solver_(linear_solver) {}

        inline Dogleg *clone() const override { return new Dogleg(std::shared_ptr<LinearSolver>(ls_solver_->clone())); }

        bool apply(const Vector &b, Vector &x) override { return aux_solve(*this->get_operator(), b, x); }

        void read(Input &in) override {
            TRSubproblem<Matrix, Vector>::read(in);

            if (ls_solver_) {
                in.get("linear-solver", *ls_solver_);
            }
        }

        void print_usage(std::ostream &os) const override {
            TRSubproblem<Matrix, Vector>::print_usage(os);
            this->print_param_usage(os, "linear-solver", "LinearSolver", "Input parameters for linear solver.", "-");
        }

        void init_memory(const Layout &layout) override {
            Bg_.zeros(layout);
            p_SD_.zeros(layout);
        }

    protected:
        bool aux_solve(const Matrix &B, const Vector &g, Vector &p_k) {
            if (ls_solver_) {
                ls_solver_->solve(B, g, p_k);
            } else {
                utopia_error("Dogleg:: linear solver not provided... \n");
            }

            if (norm2(p_k) <= this->current_radius()) {
                return true;
            } else {
                Scalar g_dots, g_B_g;
                Bg_ = B * g;

                dots(g, g, g_dots, g, Bg_, g_B_g);

                g_B_g = g_dots / g_B_g;
                p_SD_ = g_B_g * g;
                Scalar SD_norm = std::abs(g_B_g) * std::sqrt(g_dots);

                if (SD_norm >= this->current_radius()) {
                    p_k = (this->current_radius() / SD_norm) * p_SD_;
                    return true;
                } else {
                    Bg_ = p_k - p_SD_;

                    Scalar a, b;
                    dots(Bg_, Bg_, a, p_SD_, Bg_, b);
                    b *= 2.0;

                    Scalar c = (SD_norm * SD_norm) - (this->current_radius() * this->current_radius());

                    Scalar tau = this->quadratic_function(a, b, c);
                    p_k = p_SD_ + tau * Bg_;

                    return true;
                }
            }
        }

    public:
        void set_linear_solver(const std::shared_ptr<LinearSolver> &ls) { ls_solver_ = ls; }

    private:
        std::shared_ptr<LinearSolver> ls_solver_;
        Vector Bg_, p_SD_;
    };
}  // namespace utopia

#endif  // UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP
#ifndef TR_SUBPROBLEM_L2NORM_HPP
#define TR_SUBPROBLEM_L2NORM_HPP

#include "utopia_IterativeSolver.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"

#include <string>

namespace utopia {

    template <class Vector>
    class TRSubproblemBase : public virtual Configurable {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        ~TRSubproblemBase() override = default;

        TRSubproblemBase() : current_radius_(1e14) {}

        void current_radius(const Scalar &radius) { current_radius_ = radius; }

        Scalar current_radius() { return current_radius_; }

        void read(Input &in) override { in.get("current_radius", current_radius_); }

        void print_usage(std::ostream &os) const override {
            this->print_param_usage(os, "current_radius", "real", "Value of trust region radius.", "1e14");
        }

    protected:
        Scalar quad_solver(const Vector &s, const Vector &p_k, const Scalar &delta, Vector &result) {
            Scalar a, b, c, x1, x2, nom, denom, tau;

            a = dot(p_k, p_k);
            b = dot(2 * s, p_k);
            c = dot(s, s) - delta * delta;

            nom = b * b - 4 * a * c;
            nom = std::sqrt(nom);
            denom = 2 * a;

            x1 = (-b + nom) / denom;
            x2 = (-b - nom) / denom;

            tau = std::max(x1, x2);

            if (tau != tau) tau = 0;

            result = s + tau * p_k;
            return tau;
        }

        Scalar quadratic_function(const Scalar &a, const Scalar &b, const Scalar &c) {
            Scalar sqrt_discriminant = std::sqrt(b * b - 4.0 * a * c);

            Scalar lower = (-b + sqrt_discriminant) / (2.0 * a);
            Scalar upper = (-b - sqrt_discriminant) / (2.0 * a);

            return std::max(upper, lower);
        }

    protected:
        Scalar current_radius_;
    };

    template <class Matrix, class Vector>
    class TRSubproblem : public virtual IterativeSolver<Matrix, Vector>, public virtual TRSubproblemBase<Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;

    public:
        TRSubproblem() = default;

        ~TRSubproblem() override = default;
        TRSubproblem *clone() const override = 0;

        void read(Input &in) override {
            IterativeSolver<Matrix, Vector>::read(in);
            TRSubproblemBase<Vector>::read(in);
        }

        void print_usage(std::ostream &os) const override {
            IterativeSolver<Matrix, Vector>::print_usage(os);
            TRSubproblemBase<Vector>::print_usage(os);
        }
    };

    template <class Vector>
    class MatrixFreeTRSubproblem : public virtual MatrixFreeLinearSolver<Vector>,
                                   public virtual TRSubproblemBase<Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;

    public:
        MatrixFreeTRSubproblem() = default;

        ~MatrixFreeTRSubproblem() override = default;
        MatrixFreeTRSubproblem *clone() const override = 0;

        void read(Input &in) override {
            MatrixFreeLinearSolver<Vector>::read(in);
            TRSubproblemBase<Vector>::read(in);
        }

        void print_usage(std::ostream &os) const override {
            MatrixFreeLinearSolver<Vector>::print_usage(os);
            TRSubproblemBase<Vector>::print_usage(os);
        }
    };

    template <class Matrix, class Vector>
    class OperatorBasedTRSubproblem : public virtual MatrixFreeTRSubproblem<Vector>,
                                      public virtual TRSubproblem<Matrix, Vector> {
    public:
        using MatrixFreeTRSubproblem<Vector>::update;
        using TRSubproblem<Matrix, Vector>::update;
        using MatrixFreeTRSubproblem<Vector>::solve;

        ~OperatorBasedTRSubproblem() override = default;

        bool solve(const Matrix &A, const Vector &b, Vector &x) override {
            update(make_ref(A));
            return solve(operator_cast<Vector>(A), b, x);
        }

        void update(const Operator<Vector> &) override {}

        void update(const std::shared_ptr<const Matrix> &op) override {
            TRSubproblem<Matrix, Vector>::update(op);
            update(operator_cast<Vector>(*op));
        }

        bool smooth(const Vector &rhs, Vector &x) override {
            SizeType temp = this->max_it();
            this->max_it(this->sweeps());
            solve(operator_cast<Vector>(*this->get_operator()), rhs, x);
            this->max_it(temp);
            return true;
        }

        bool apply(const Vector &b, Vector &x) override {
            return solve(operator_cast<Vector>(*this->get_operator()), b, x);
        }

        OperatorBasedTRSubproblem *clone() const override = 0;

        void read(Input &in) override {
            MatrixFreeTRSubproblem<Vector>::read(in);
            TRSubproblem<Matrix, Vector>::read(in);
        }

        void print_usage(std::ostream &os) const override {
            MatrixFreeTRSubproblem<Vector>::print_usage(os);
            TRSubproblem<Matrix, Vector>::print_usage(os);
        }
    };

}  // namespace utopia

#endif  // TR_SUBPROBLEM_L2NORM_HPP

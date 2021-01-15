#ifndef UTOPIA_MORE_SORENSEN_EIGEN_QP_SUBPROBLEM_HPP
#define UTOPIA_MORE_SORENSEN_EIGEN_QP_SUBPROBLEM_HPP

#include "utopia_EigenSolver.hpp"
#include "utopia_Layout.hpp"
#include "utopia_TRSubproblem.hpp"

namespace utopia {
    /**
     * @brief      Class for More Sorensen minimization algorithm, where
     * initialization of lambda_0 is based on eigen sol.
     *
     * 			   WARNING:: 	Computation of direction could be more
     * efficient by using cholesky decomposition and store just factors, but it is a
     * bit anoying to do so with petsc => we solve system 2 \times, which is less
     * efficient, but ok as this is just proof of concept solver
     */
    template <class Matrix, class Vector>
    class MoreSorensenEigen final : public TRSubproblem<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;
        using EigenSolver = utopia::EigenSolver<Matrix, Vector>;

    public:
        MoreSorensenEigen(const std::shared_ptr<LinearSolver> &linear_solver,
                          const std::shared_ptr<EigenSolver> &eigen_solver)
            : TRSubproblem<Matrix, Vector>(),
              linear_solver_(linear_solver),
              eigen_solver_(eigen_solver),
              kappa_easy_(1e-10),
              lambda_eps_(1e-5),
              initialized_(false) {}

        inline void kappa_easy(const Scalar &kappa) { kappa_easy_ = kappa; }

        inline Scalar kappa_easy() { return kappa_easy_; }

        inline void lambda_eps(const Scalar &lambda_eps) { lambda_eps_ = lambda_eps; }

        inline Scalar lambda_eps() { return lambda_eps_; }

        virtual void set_linear_solver(const std::shared_ptr<LinearSolver> &ls) { linear_solver_ = ls; }

        MoreSorensenEigen *clone() const override { return new MoreSorensenEigen(*this); }

        bool apply(const Vector &b, Vector &x) override {
            auto b_layout = layout(b);
            if (!initialized_ || !b_layout.same(layout_)) {
                init_memory(b_layout);
            }

            return aux_solve(*this->get_operator(), b, x);
        }

        void read(Input &in) override {
            TRSubproblem<Matrix, Vector>::read(in);
            in.get("kappa_easy", kappa_easy_);
            in.get("lambda_eps", lambda_eps_);

            if (linear_solver_) {
                in.get("linear-solver", *linear_solver_);
            }
            if (eigen_solver_) {
                in.get("eigen-solver", *eigen_solver_);
            }
        }

        void print_usage(std::ostream &os) const override {
            TRSubproblem<Matrix, Vector>::print_usage(os);

            this->print_param_usage(os, "kappa_easy", "double", "Treshold for kappa estimation.", "1e-10");
            this->print_param_usage(os, "lambda_eps", "double", "Treshold for min. eigenvalue.", "1e-5");

            this->print_param_usage(os, "linear-solver", "LinearSolver", "Input parameters for linear solver.", "-");
            this->print_param_usage(os, "eigen-solver", "EigenSolver", "Input parameters for eigen solver.", "-");
        }

    private:
        bool aux_solve(const Matrix &H, const Vector &g, Vector &s_k) {
            Scalar lambda, s_norm, lambda_old;

            // UTOPIA_NO_ALLOC_BEGIN("MoreSorensen:region1");
            if (empty(s_k) || size(s_k) != size(g))
                s_k = 0.0 * g;
            else
                s_k.set(0.0);

            // ---------------------- initialization  of lambda_0
            // ------------------------
            eigen_solver_->portion_of_spectrum("smallest_real");
            eigen_solver_->number_of_eigenvalues(1);
            eigen_solver_->solve(H);

            eigen_solver_->get_real_eigenpair(0, lambda, eigenvector_);

            // 	decide if PD case or not
            lambda = (lambda > 0.0) ? 0.0 : -1.0 * (lambda - lambda_eps_);

            H_lambda_ = H;
            if (lambda != 0.0) {
                H_lambda_.shift_diag(lambda);
            }

            linear_solver_->solve(H_lambda_, g, s_k);
            s_norm = norm2(s_k);
            // UTOPIA_NO_ALLOC_END();

            if (s_norm <= this->current_radius()) {
                if (lambda == 0.0 || s_norm == this->current_radius())
                    return true;
                else {
                    // UTOPIA_NO_ALLOC_BEGIN("MoreSorensen:region2");
                    Scalar s_k_eigen, sk_sk;
                    dots(s_k, eigenvector_, s_k_eigen, s_k, s_k, sk_sk);

                    // we are in hard case, let's find solution on boundary, which is
                    // orthogonal to E_1
                    //                     because eigenvector is normalized
                    Scalar alpha = this->quadratic_function(
                        1.0, 2.0 * s_k_eigen, sk_sk - (this->current_radius() * this->current_radius()));
                    s_k += alpha * eigenvector_;
                    // UTOPIA_NO_ALLOC_END();
                    return true;
                }
            }

            for (auto it = 0; it < this->max_it(); it++) {
                if (std::abs(s_norm - this->current_radius()) <= kappa_easy_ * this->current_radius()) return true;

                lambda_old = lambda;

                // UTOPIA_NO_ALLOC_BEGIN("MoreSorensen:region3");
                grad_s_lambda_.set(0);

                s_k *= -1.0;
                linear_solver_->solve(H_lambda_, s_k, grad_s_lambda_);

                Scalar grad = 1.0 / s_norm - 1.0 / this->current_radius();
                Scalar hessian = dot(s_k, grad_s_lambda_) / std::pow(s_norm, 3);

                lambda -= grad / hessian;

                // H should be H + \lambda I
                H_lambda_.shift_diag(lambda - lambda_old);

                s_k.set(0.0);
                linear_solver_->solve(H_lambda_, g, s_k);
                s_norm = norm2(s_k);
                // UTOPIA_NO_ALLOC_END();
            }

            return true;
        }

        void init_memory(const Layout &layout) override {
            // resets all buffers in case the size has changed
            eigenvector_.zeros(layout);
            grad_s_lambda_.zeros(layout);

            initialized_ = true;
            layout_ = layout;
        }

    private:
        std::shared_ptr<LinearSolver> linear_solver_;
        std::shared_ptr<EigenSolver> eigen_solver_;

        Scalar kappa_easy_;
        Scalar lambda_eps_;

        bool initialized_;
        Layout layout_;

        Vector eigenvector_;
        Vector grad_s_lambda_;

        Matrix H_lambda_;
    };

}  // namespace utopia

#endif  // UTOPIA_MORE_SORENSEN_EIGEN_QP_SUBPROBLEM_HPP
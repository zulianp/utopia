#ifndef UTOPIA_BRATU_1D_HPP
#define UTOPIA_BRATU_1D_HPP

#include "utopia.hpp"
#include "utopia_TestFunctions.hpp"

namespace utopia {

    template <typename Matrix, typename Vector>
    class Bratu1D final : virtual public UnconstrainedExtendedTestFunction<Matrix, Vector>,
                          virtual public ConstrainedExtendedTestFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        Bratu1D(const SizeType &n) : n_(n), lambda_(3.5), lambda_critical_(3.513830719) {
            init(Comm::get_default(), n);
        }

        Bratu1D(const Comm &comm = Comm::get_default(), const SizeType &n = 10)
            : n_(n), lambda_(3.5), lambda_critical_(3.513830719) {
            init(comm, n);
        }

        void init(const Comm &comm, const SizeType &n) {
            n_ = n;
            x0_.zeros(layout(comm, Traits::decide(), n));
            auto vec_layout = layout(x0_);
            exact_sol_.values(vec_layout, 0.0);
            A_help_ = make_unique<Vector>(vec_layout, 0.0);

            a_ = 0.0;
            b_ = 1.0;
            L_ = b_ - a_;
            h_ = L_ / (n_ - 1);

            // check if user param reasonable
            if (lambda_ > lambda_critical_) lambda_ = 3.5;

            H_.sparse(square_matrix_layout(vec_layout), 3, 2);
            assemble_laplacian_1D(H_);

            // bc_indices_.push_back(0);
            // bc_indices_.push_back(n-1);

            Vector bc_markers(vec_layout, 0.0);
            {
                Write<Vector> wv(bc_markers);

                Range r = range(bc_markers);

                if (r.begin() == 0) {
                    bc_markers.set(0, 1.0);
                    bc_indices_.push_back(0);
                }

                if (r.end() == n_) {
                    bc_markers.set(n - 1, 1.0);
                    bc_indices_.push_back(n - 1);
                }
            }

            // disp(bc_markers);
            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, x0_);

            // this->constraints_ = make_box_constaints(std::make_shared<Vector>(values(n_, -9e9)),
            //                                          std::make_shared<Vector>(values(n_, 9e9)));
        }

        bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const override {
            H = H_;
            set_zero_rows(H, bc_indices_, 1.);
            return true;
        }

        bool value(const Vector &x, Scalar &energy) const override {
            *A_help_ = ((H_)*x);
            energy = 0.5 * dot(x, *A_help_);
            *A_help_ = exp(x);
            energy -= h_ * lambda_ * sum(*A_help_);

            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override {
            // UTOPIA_NO_ALLOC_BEGIN("Bratu1D::gradient1");
            g = (H_ * x);
            // UTOPIA_NO_ALLOC_END();

            // UTOPIA_NO_ALLOC_BEGIN("Bratu1D::gradient2");
            *A_help_ = h_ * lambda_ * exp(x);
            g -= *A_help_;
            // UTOPIA_NO_ALLOC_END();

            // UTOPIA_NO_ALLOC_BEGIN("Bratu1D::gradient3");
            {
                Write<Vector> w(g);
                Read<Vector> read(x);

                Range r = range(g);

                if (r.begin() == 0) {
                    g.set(0, -x.get(0));
                }
                if (r.end() == n_) {
                    g.set(n_ - 1, -x.get(n_ - 1));
                }
            }
            // UTOPIA_NO_ALLOC_END();

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override {
            H = H_;
            *A_help_ = h_ * (-lambda_) * exp(x);
            H += Matrix(diag(*A_help_));
            set_zero_rows(H, bc_indices_, 1.);

            return true;
        }

        bool hessian(const Vector &x, Matrix &result, Matrix &prec) const override {
            UTOPIA_UNUSED(x);
            UTOPIA_UNUSED(result);
            UTOPIA_UNUSED(prec);
            return false;
        }

        bool has_preconditioner() const override { return false; }

        Vector initial_guess() const override { return x0_; }

        const Vector &exact_sol() const override { return exact_sol_; }

        Scalar min_function_value() const override {
            // depends on the solution to which we converged to
            std::cout << "Bratu1D:: min_function_value :: wrong.... \n";
            return -1.012;
        }

        std::string name() const override { return "Bratu1D"; }

        SizeType dim() const override { return n_; }

        bool exact_sol_known() const override { return false; }

        bool parallel() const override { return true; }

        void apply_bc_to_initial_guess(Vector &x) {
            // here, we assume that BC are 0 on the initial and end part of subdomain
            {
                Write<Vector> w(x);
                Range r = range(x);

                if (r.begin() == 0) {
                    x.set(0, 0.0);
                }

                if (r.end() == n_) {
                    x.set(n_ - 1, 0.0);
                }
            }
        }

        void generate_constraints(Vector &lb,
                                  Vector &ub,
                                  const Scalar lower_const = -1.0 * std::numeric_limits<Scalar>::infinity(),
                                  const Scalar upper_const = std::numeric_limits<Scalar>::infinity()) {
            lb.values(layout(x0_), lower_const);
            ub.values(layout(x0_), upper_const);
        }

    private:
        void assemble_laplacian_1D(Matrix &M) {
            {
                // n x n matrix with maximum 3 entries x row
                Write<Matrix> w(M);
                Range r = row_range(M);
                auto n = size(M).get(0);

                for (SizeType i = r.begin(); i != r.end(); ++i) {
                    if (i > 0) {
                        M.set(i, i - 1, -1.0);
                    }

                    if (i < n - 1) {
                        M.set(i, i + 1, -1.0);
                    }

                    if (i == 0 || i == n - 1) {
                        M.set(i, i, 1.);
                    } else {
                        M.set(i, i, 2.0);
                    }
                }
            }

            // auto n = size(M).get(0);
            // M *= 1./(h_*h_);
            M *= 1. / (h_);
        }

    private:
        Scalar n_, L_, h_;
        Scalar lambda_;
        const Scalar lambda_critical_;
        Scalar a_, b_;

        std::vector<SizeType> bc_indices_;

        Matrix H_;
        Vector x0_;
        Vector exact_sol_;

        std::unique_ptr<Vector> A_help_;
    };

}  // namespace utopia

#endif  // UTOPIA_BRATU_1D_HPP
#ifndef UTOPIA_POWER_METHOD_HPP
#define UTOPIA_POWER_METHOD_HPP

#include <memory>
#include "utopia_MatrixFreeLinearSolver.hpp"

namespace utopia {

    // TODO:: see if there is more appropriate interface to be used...
    /**
     * @brief      Power method
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Vector>
    class PowerMethod final : public virtual Clonable, public Configurable {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

    public:
        PowerMethod() : tol_(1e-5), max_it_(5) {}

        PowerMethod(const PowerMethod &other) : tol_(other.tol_), max_it_(other.max_it_) {}

        void read(Input &in) override { in.get("use_rand_vec_init", use_rand_vec_init_); }

        void tol(const Scalar &tol) { tol_ = tol; }
        Scalar tol() const { return tol_; }

        void max_it(const SizeType &power_method_max_it) { max_it_ = power_method_max_it; }
        SizeType max_it() const { return max_it_; }

        void use_rand_vec_init(const bool &use_rand_vec_init) { use_rand_vec_init_ = use_rand_vec_init; }
        bool use_rand_vec_init() const { return use_rand_vec_init_; }

        void init_memory(const Layout &layout) {
            assert(layout.local_size() > 0);

            // resets all buffers in case the size has changed
            help_f1.zeros(layout);
            eigenvector_.zeros(layout);
            fi_.zeros(layout);

            initialized_ = true;
            layout_ = layout;
        }

        void print_usage(std::ostream &os) const override {
            // OperatorBasedLinearSolver<Matrix, Vector>::print_usage(os);
        }

        PowerMethod *clone() const override { return new PowerMethod(*this); }

        // void copy(const PowerMethod &other) { Super::operator=(other); }

        Scalar get_max_eig(const Operator<Vector> &A, const Vector &rhs, const Scalar &shift = 0.0) {
            // Super simple power method to estimate the biggest eigenvalue
            assert(!empty(eigenvector_));

            if (use_rand_vec_init_ == true) {
                auto d_eigenvector_ = local_view_device(eigenvector_);

                parallel_for(
                    local_range_device(eigenvector_), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar val = ((Scalar)std::rand() / (RAND_MAX)) + 1;

                        d_eigenvector_.set(i, val);
                    });
            } else if (norm2(eigenvector_) == 0) {
                eigenvector_ = rhs;
            }

            get_max_eig(A, shift);
        }

        Scalar get_max_eig(const Operator<Vector> &A, const Scalar &shift = 0.0) {
            // normalize IG
            eigenvector_ = Scalar(1. / norm2(eigenvector_)) * eigenvector_;

            SizeType it = 0;
            bool converged = false;
            Scalar gnorm, lambda = 0.0, lambda_old;

            while (!converged) {
                help_f1 = eigenvector_;
                A.apply(help_f1, eigenvector_);

                // todo:: verify
                if (shift != 0.0) {
                    eigenvector_ - shift *help_f1;
                }

                eigenvector_ = Scalar(1.0 / Scalar(norm2(eigenvector_))) * eigenvector_;

                lambda_old = lambda;

                A.apply(eigenvector_, help_f1);
                lambda = dot(eigenvector_, help_f1);

                fi_ = eigenvector_ - help_f1;
                gnorm = norm2(fi_);

                converged = ((gnorm < tol_) || (std::abs(lambda_old - lambda) < tol_) || it > max_it_) ? true : false;

                it = it + 1;
            }

            if (this->verbose() && mpi_world_rank() == 0)
                utopia::out() << "Power method converged in " << it << " iterations. Largest eig: " << lambda << "  \n";

            return lambda;
        }

        bool initialized_{false};
        Layout layout_;
        Scalar tol_;
        SizeType max_it_;

        bool use_rand_vec_init_{false};

        // This fields are not to be copied anywhere
        Vector help_f1, eigenvector_, fi_;
    };
}  // namespace utopia

#endif

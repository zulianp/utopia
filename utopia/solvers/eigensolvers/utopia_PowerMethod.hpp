#ifndef UTOPIA_POWER_METHOD_HPP
#define UTOPIA_POWER_METHOD_HPP

#include <memory>
#include "utopia_MatrixFreeLinearSolver.hpp"

namespace utopia {

    // TODO:: see if there is more appropriate interface to be used... such as EigenSolver...
    /**
     * @brief      Power method
     * @tparam     Vector
     */
    template <class Vector>
    class PowerMethod final : public virtual Clonable, public Configurable {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

    public:
        PowerMethod() {}

        PowerMethod(const PowerMethod &other) : tol_(other.tol_), max_it_(other.max_it_) {}

        void read(Input &in) override {
            in.get("use_rand_vec_init", use_rand_vec_init_);
            in.get("verbose", verbose_);
            in.get("tol", tol_);
            in.get("max_it", max_it_);
        }

        bool verbose() const { return verbose_; }
        void verbose(const bool &flg) { verbose_ = flg; }

        void tol(const Scalar &tol) { tol_ = tol; }
        Scalar tol() const { return tol_; }

        void max_it(const SizeType &power_method_max_it) { max_it_ = power_method_max_it; }
        SizeType max_it() const { return max_it_; }

        void use_rand_vec_init(const bool &use_rand_vec_init) { use_rand_vec_init_ = use_rand_vec_init; }
        bool use_rand_vec_init() const { return use_rand_vec_init_; }

        void init_memory(const Layout &layout) {
            assert(layout.local_size() > 0);

            // resets all buffers in case the size has changed
            eigenvector_.zeros(layout);

            initialized_ = true;
            layout_ = layout;
        }

        void print_usage(std::ostream & /*os*/) const override {}

        PowerMethod *clone() const override { return new PowerMethod(*this); }

        Scalar get_max_eig(const Operator<Vector> &A, const Vector &rhs) {
            // Super simple power method to estimate the largest eigenvalue
            assert(!empty(eigenvector_));

            if (use_rand_vec_init_ == true) {
                auto d_eigenvector_ = local_view_device(eigenvector_);

                parallel_for(
                    local_range_device(eigenvector_), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar val = ((Scalar)std::rand() / (RAND_MAX)) + 1;

                        d_eigenvector_.set(i, val);
                    });
            }
            // for first it, we use IG for eigenvector to be rhs
            // for all next its, we use eigenvector computed on previous run
            // unless specified to be random, then we reinit
            else if (norm2(eigenvector_) == 0) {
                eigenvector_ = rhs;
            }

            return compute_max_eig(A);
        }

        Scalar get_max_eig(const Operator<Vector> &A) {
            // Super simple power method to estimate the largest eigenvalue
            assert(!empty(eigenvector_));

            if (use_rand_vec_init_ == true || norm2(eigenvector_) == 0) {
                auto d_eigenvector_ = local_view_device(eigenvector_);

                parallel_for(
                    local_range_device(eigenvector_), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar val = ((Scalar)std::rand() / (RAND_MAX)) + 1;

                        d_eigenvector_.set(i, val);
                    });
            }

            return compute_max_eig(A);
        }

    private:
        Scalar compute_max_eig(const Operator<Vector> &A) {
            // normalize IG
            eigenvector_ = Scalar(1. / norm2(eigenvector_)) * eigenvector_;

            SizeType it = 0;
            bool converged = false;
            Scalar lambda = 0.0, lambda_old = 0.0;

            while (!converged) {
                A.apply(eigenvector_, eigenvector_);

                lambda = norm2(eigenvector_);
                eigenvector_ = (1. / lambda) * eigenvector_;

                converged = ((device::abs(lambda_old - lambda) < tol_) || it > max_it_) ? true : false;
                lambda_old = lambda;

                it = it + 1;
            }

            if (this->verbose() && mpi_world_rank() == 0) {
                utopia::out() << "Power method converged in " << it << " iterations. Largest eig: " << lambda << "  \n";
            }
            return lambda;
        }

        bool initialized_{false};
        Layout layout_;
        Scalar tol_{1e-2};
        SizeType max_it_{5};
        bool verbose_{false};
        bool use_rand_vec_init_{false};

        // This fields are not to be copied anywhere
        Vector eigenvector_;
    };
}  // namespace utopia

#endif

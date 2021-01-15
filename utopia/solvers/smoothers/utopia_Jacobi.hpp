#ifndef UTOPIA_JACOBI_HPP
#define UTOPIA_JACOBI_HPP
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"

namespace utopia {

    /**
     * @brief Good example on how to implement a Linear Solver, Preconditioner, and
     * Smoother within utopia. Precomputations are done in the update method, and,
     * in order to avoid costly allocations, possible temporaries are stored as
     * member variables.
     */
    template <class Matrix, class Vector>
    class Jacobi final : public IterativeSolver<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        typedef utopia::IterativeSolver<Matrix, Vector> Solver;

    public:
        /**
         * @brief      Very simple Jacobi solver.
         *
         * @param[in]  omega  The relaxation parameter (unused atm).
         */
        Jacobi() : compute_norm_each_(10), preconditioner_mode_(false) {}

        void read(Input &in) override { Solver::read(in); }

        void print_usage(std::ostream &os) const override { Solver::print_usage(os); }

        void compute_norm_each(const SizeType &n_iter) { compute_norm_each_ = n_iter; }

        void preconditioner_mode(const bool &val) { preconditioner_mode_ = val; }

        bool apply(const Vector &rhs, Vector &x) override {
            const Matrix &A = *this->get_operator();

            SizeType it = 0;
            UTOPIA_NO_ALLOC_BEGIN("Jacobi:r2");
            r_ = rhs - A * x;
            UTOPIA_NO_ALLOC_END();

            Scalar g_norm0 = 0.0;

            if (!preconditioner_mode_) {
                g_norm0 = norm2(r_);
            }

            Scalar g_norm = g_norm0;

            this->init_solver("Point Jacobi", {"it. ", "||r||"});

            // First step for using r0
            x += e_mul(d_inv_, r_);

            while (true) {
                it++;

                if (preconditioner_mode_) {
                    if (it >= this->max_it()) {
                        // g_norms are just zero
                        return this->check_convergence(it, g_norm, g_norm / g_norm0, 1);
                    }
                } else {
                    if (it % compute_norm_each_ == 0 || it >= this->max_it()) {
                        g_norm = norm2(r_);

                        if (this->verbose()) {
                            PrintInfo::print_iter_status(it, {g_norm});
                        }

                        if (this->check_convergence(it, g_norm, g_norm / g_norm0, 1)) {
                            return true;
                        }
                    }
                }

                r_ = rhs - A * x;
                x += e_mul(d_inv_, r_);
            }

            return false;
        }

        bool smooth(const Vector &rhs, Vector &x) override {
            for (SizeType it = 0; it < this->sweeps(); it++) {
                sweep(rhs, x);
            }

            return true;
        }

        inline Jacobi *clone() const override { return new Jacobi(*this); }

        void update(const std::shared_ptr<const Matrix> &op) override {
            Solver::update(op);
            const auto &A = *op;

            UTOPIA_NO_ALLOC_BEGIN("Jacobi:r4");
            d_inv_ = diag(A);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("Jacobi:r4.2");
            e_pseudo_inv(d_inv_, d_inv_, 0);
            UTOPIA_NO_ALLOC_END();
        }

        void init_memory(const Layout &layout) override {
            d_inv_.zeros(layout);
            r_.zeros(layout);
        }

    private:
        Vector d_inv_;
        Vector r_;
        SizeType compute_norm_each_;
        bool preconditioner_mode_;

        /**
         * @brief      Checks if there is a zero in the vector, if yes turn it into 1.
         *
         * @param      diag_A  { D_{-1}}
         *
         * @return     {  }
         */

        inline bool sweep(const Vector &rhs, Vector &x) {
            UTOPIA_NO_ALLOC_BEGIN("Jacobi:r1");

            const Matrix &A = *this->get_operator();
            r_ = rhs - A * x;
            x += e_mul(d_inv_, r_);

            UTOPIA_NO_ALLOC_END();
            return true;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_JACOBI_HPP

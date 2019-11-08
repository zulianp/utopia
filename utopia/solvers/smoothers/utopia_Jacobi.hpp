#ifndef UTOPIA_JACOBI_HPP
#define UTOPIA_JACOBI_HPP
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"

namespace utopia {

    /**
     * @brief Good example on how to implement a Linear Solver, Preconditioner, and Smoother
     * within utopia. Precomputations are done in the update method, and, in order to avoid
     * costly allocations, possible temporaries are stored as member variables.
     */
    template<class Matrix, class Vector>
    class Jacobi final: public Smoother<Matrix, Vector>, public IterativeSolver<Matrix, Vector> {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::IterativeSolver<Matrix, Vector> Solver;
        typedef utopia::Smoother<Matrix, Vector> Smoother;

    public:
        /**
         * @brief      Very simple Jacobi solver.
         *
         * @param[in]  omega  The relaxation parameter (unused atm).
         */
        Jacobi()
        {

        }

        void read(Input &in) override
        {
            Solver::read(in);
            Smoother::read(in);
        }

        void print_usage(std::ostream &os) const override
        {
            Solver::print_usage(os);
            Smoother::print_usage(os);
        }


        bool apply(const Vector &rhs, Vector &x) override
        {
            const Matrix &A = *this->get_operator();

            SizeType it = 0;
            UTOPIA_NO_ALLOC_BEGIN("Jacobi:r2");
            r_ = rhs - A * x;
            UTOPIA_NO_ALLOC_END();

            Scalar g_norm0 = norm2(r_);
            Scalar g_norm = g_norm0;
            SizeType compute_norm_each = 50;

            this->init_solver("Point Jacobi", {"it. ", "||r||" });

            //First step for using r0
            x += e_mul(d_inv_, r_);

            while(true) {
                
                if(it++ % compute_norm_each == 0) {
                    g_norm = norm2(r_);

                    if(this->verbose()) {
                        PrintInfo::print_iter_status(it, {g_norm});
                    }

                    if(this->check_convergence(it, g_norm, g_norm/g_norm0, 1)) {
                        return true;
                    }
                }

                r_ = rhs - A * x;
                x += e_mul(d_inv_, r_);
            }

            return false;
        }

        bool smooth(const Vector &rhs, Vector &x) override
        {
            for(SizeType it = 0; it < this->sweeps(); it++)
            {
                sweep(rhs, x);
            }

            return true;
        }

        inline Jacobi * clone() const override
        {
            return new Jacobi(*this);
        }

        void update(const std::shared_ptr<const Matrix> &op) override
        {
            Solver::update(op);
            const auto &A = *op;

            UTOPIA_NO_ALLOC_BEGIN("Jacobi:r4");
            d_inv_ = diag(A);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("Jacobi:r4.2");
            e_pseudo_inv(d_inv_, d_inv_, 0);
            UTOPIA_NO_ALLOC_END();
        }

    private:
        Vector d_inv_;
        Vector r_;

        /**
         * @brief      Checks if there is a zero in the vector, if yes turn it into 1.
         *
         * @param      diag_A  { D_{-1}}
         *
         * @return     {  }
         */

        inline bool sweep(const Vector &rhs, Vector &x)
        {
            UTOPIA_NO_ALLOC_BEGIN("Jacobi:r1");

            const Matrix &A = *this->get_operator();
            r_ = rhs - A * x;
            x += e_mul(d_inv_, r_);

            UTOPIA_NO_ALLOC_END();
            return true;
        }

    };

}

#endif //UTOPIA_JACOBI_HPP


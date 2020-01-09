#ifndef UTOPIA_POINT_JACOBI_HPP
#define UTOPIA_POINT_JACOBI_HPP
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"

namespace utopia {

    /**
     * @brief Good example on how to implement a Linear Solver, Preconditioner, and Smoother
     * within utopia. Precomputations are done in the update method, and, in order to avoid
     * costly allocations, possible temporaries are stored as member variables.
     */
    template<class Matrix, class Vector>
    class PointJacobi final: public IterativeSolver<Matrix, Vector> 
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::IterativeSolver<Matrix, Vector> Solver;

    public:
        /**
         * @brief      Very simple Jacobi solver.
         *
         * @param[in]  omega  The relaxation parameter (unused atm).
         */
        PointJacobi()
        {

        }

        void read(Input &in) override
        {
            Solver::read(in);
        }

        void print_usage(std::ostream &os) const override
        {
            Solver::print_usage(os);
        }


        bool apply(const Vector &rhs, Vector &x) override
        {
            const Matrix &A = *this->get_operator();

            SizeType it = 0;
            UTOPIA_NO_ALLOC_BEGIN("PointJacobi:r2");
            r_ = rhs - A * x;
            UTOPIA_NO_ALLOC_END();

            Scalar g_norm0 = norm2(r_);
            Scalar g_norm = g_norm0;
            SizeType compute_norm_each = 50;

            this->init_solver("Point Jacobi", {"it. ", "||r||" });

            while(true) {
                sweep(rhs, x);

                if(it++ % compute_norm_each == 0) {
                    UTOPIA_NO_ALLOC_BEGIN("PointJacobi:r21");
                    r_ = rhs - A * x;
                    g_norm = norm2(r_);
                    UTOPIA_NO_ALLOC_END();

                    if(this->verbose()) {
                        PrintInfo::print_iter_status(it, {g_norm});
                    }

                    if(this->check_convergence(it, g_norm, g_norm/g_norm0, 1)) {
                        return true;
                    }
                }
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

        inline PointJacobi * clone() const override
        {
            return new PointJacobi(*this);
        }

        void update(const std::shared_ptr<const Matrix> &op) override
        {
            Solver::update(op);

            const auto &A = *op;

            UTOPIA_NO_ALLOC_BEGIN("PointJacobi:r4");
            d_inv_ = diag(A);
            UTOPIA_NO_ALLOC_END();

            
            // lower and upper part of A
            LU_ = A;
            UTOPIA_NO_ALLOC_BEGIN("PointJacobi:r4.1");
            LU_ -= Matrix(diag(d_inv_));
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("PointJacobi:r4.2");
            d_inv_ = 1. / d_inv_;
            UTOPIA_NO_ALLOC_END();

            // prevents system from being indefinite
            check_indef(d_inv_);

        }


        void init_memory(const SizeType & ls) override
        {
            auto zero_expr = local_zeros(ls);
            r_ = zero_expr;
            d_inv_ = zero_expr;            
        }

    private:
        Vector d_inv_;
        Matrix LU_;
        Vector r_;

        /**
         * @brief      Checks if there is a zero in the vector, if yes turn it into 1.
         *
         * @param      diag_A  { D_{-1}}
         *
         * @return     {  }
         */
        inline bool check_indef(Vector &diag_A)
        {
            ReadAndWrite<Vector> w(diag_A);
            Range rr = range(diag_A);

            for (SizeType i = rr.begin(); i != rr.end(); i++)
            {
                if(diag_A.get(i) == 0)
                {
                    diag_A.set(i, 1);
                }
            }

            return true;
        }


        inline bool sweep(const Vector &rhs, Vector &x)
        {
            UTOPIA_NO_ALLOC_BEGIN("PointJacobi:r1");
            r_ = rhs - (LU_ * x);
            x = e_mul(d_inv_, r_);

            // const Matrix &A = *this->get_operator();
            // r_ = rhs - A * x;
            // x += e_mul(d_inv_, r_);

            UTOPIA_NO_ALLOC_END();
            return true;
        }

    };

}

#endif //UTOPIA_POINT_JACOBI_HPP


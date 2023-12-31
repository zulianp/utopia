#ifndef UTOPIA_PETSC_GS_HPP
#define UTOPIA_PETSC_GS_HPP

#include "utopia_Core.hpp"
#include "utopia_GaussSeidel.hpp"
#include "utopia_Smoother.hpp"

#include <petscksp.h>
#include <petscpc.h>

// extern "C"
// {
#include "petscmat.h"
#include "petscvec.h"
// }

namespace utopia {

    /**
     * @brief      Wrapper for PETSC implementation of SOR.
     *             Be aware, that this function doesn't run in parallel.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class GaussSeidel<Matrix, Vector, PETSC> : public IterativeSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        typedef utopia::IterativeSolver<Matrix, Vector> Solver;

    public:
        GaussSeidel() : n_local_sweeps_(5) {}

        void read(Input &in) override {
            Solver::read(in);
            in.get("n_local_sweeps", n_local_sweeps_);
        }

        void print_usage(std::ostream &os) const override { Solver::print_usage(os); }

        /**
         * @brief      Smoothing of GS from Petsc. Currently we are using symmetric block GS (builds block jacobi and on
         * blocks calls GS).
         *
         * @param[in]  A     The stiffness matrix.
         * @param[in]  rhs   The right hand side.
         * @param      x     The solution.
         */
        bool smooth(const Vector &rhs, Vector &x) override {
            const Matrix &A = *this->get_operator();

            if (this->verbose()) {
                utopia::out() << "before smooth: " << double(norm2(rhs - A * x)) << std::endl;
            }

            MatSOR(raw_type(A),
                   raw_type(rhs),
                   1,
                   //  SOR_FORWARD_SWEEP,
                   SOR_LOCAL_SYMMETRIC_SWEEP,
                   0,
                   this->sweeps(),
                   n_local_sweeps_,
                   raw_type(x));

            if (this->verbose()) {
                utopia::out() << "after smooth: " << double(norm2(rhs - A * x)) << std::endl;
            }

            return true;
        }

        /**
         * @brief      Solving system with Gauss-Seidel method.
         *
         * @param[in]  rhs   The right hand side.
         * @param      x     The solution.
         */
        bool apply(const Vector &rhs, Vector &x) override {
            const Matrix &A = *this->get_operator();

            SizeType it = 0;
            Scalar r_norm = 9999;
            this->init_solver("Petsc Gauss-Seidel", {"it. ", "||r||"});

            while (it < this->max_it() && r_norm > this->rtol()) {
                MatSOR(raw_type(A),
                       raw_type(rhs),
                       1,
                       // SOR_FORWARD_SWEEP,
                       SOR_LOCAL_SYMMETRIC_SWEEP,  // parallel implementation - builds block jacobi and on blocks it
                                                   // calls GS
                       0,
                       this->sweeps(),
                       n_local_sweeps_,
                       raw_type(x));

                r_norm = norm2(A * x - rhs);

                it += this->sweeps();

                if (this->verbose()) PrintInfo::print_iter_status(it, {r_norm});
            }

            return true;
        }

        inline GaussSeidel *clone() const override { return new GaussSeidel(*this); }

        void update(const std::shared_ptr<const Matrix> &op) override { Solver::update(op); }

    private:
        SizeType n_local_sweeps_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_GS_HPP

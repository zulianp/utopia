#ifndef UTOPIA_SP_STATIC_CONDENSATION_KRYLOV_HPP
#define UTOPIA_SP_STATIC_CONDENSATION_KRYLOV_HPP

#include <memory>
#include "utopia_LinearSolver.hpp"
#include "utopia_BiCGStab.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"

namespace utopia {

    /**
    @brief Solves the following saddle-point problem:
        | A_m	+ T^T A_s T |  u_m = | f_m + T^Tf_s |
         u_s  = | T u_m |
    m := master, s := slave
    **/

    template<class Matrix, class Vector>
    class SPStaticCondensationKrylov {
    public:

        SPStaticCondensationKrylov()
        : linear_solver_(std::make_shared<BiCGStab<Matrix, Vector, HOMEMADE>>()), coupling_op_solver_(std::make_shared<BiCGStab<Matrix, Vector>>())
        {}

        void coupling_op_solver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &solver)
        {
            coupling_op_solver_ = solver;
        }

        //h
        void linear_solver(const std::shared_ptr<MatrixFreeLinearSolver<Vector>> &solver)
        {
            linear_solver_ = solver;
        }

        void update(
            const std::shared_ptr<Matrix> &A_m,
            const std::shared_ptr<Matrix> &A_s,
            const std::shared_ptr<Matrix> &B,
            const std::shared_ptr<Matrix> &D)
        {
            set_up(B, D);
            update(A_m, A_s);
        }

        void set_up(
            const std::shared_ptr<Matrix> &B,
            const std::shared_ptr<Matrix> &D)
        {
            this->B = B;
            this->D = D;

            coupling_op_solver_->update(D);
        }

        void update(
            const std::shared_ptr<Matrix> &A_m,
            const std::shared_ptr<Matrix> &A_s)
        {
            this->A_m = A_m;
            this->A_s = A_s;

            // S = utopia::op<Vector>(
            //     A_m->comm(),
            //     A_m->size(),
            //     A_m->local_size(),
            //     [this](const Vector &b, Vector &x) -> bool
            //     {
            //         apply_T(b, x);
            //         Vector x_temp = (*this->A_s) * x;
            //         apply_T_transpose(x_temp, x);
            //         x_temp = (*this->A_m) * b;
            //         x += x_temp;
            //         return true;
            //     });
        }

        bool apply(
            const Vector &rhs_m,
            const Vector &rhs_s,
            Vector &sol_m,
            Vector &sol_s
        )
        {
            assert(!empty(rhs_m));
            assert(!empty(rhs_s));

            if(!empty(rhs)) {
                rhs.set(0.);
            }

            if(empty(sol_m)) {
                sol_m = local_zeros(local_size(rhs_m));
            }

            if(empty(sol_s)) {
                sol_s = local_zeros(local_size(rhs_s));
            }

            if(!apply_T_transpose(rhs_s, rhs)) return false;
            rhs += rhs_m;

            bool ok = linear_solver_->solve(*S, rhs, sol_m);
            //even if failed create the slave solution vector
            return apply_T(sol_m, sol_s) && ok;
        }

    private:
        Vector rhs;
        std::shared_ptr<Matrix> A_m, A_s, B, D;
        std::shared_ptr<MatrixFreeLinearSolver<Vector>> linear_solver_;
        std::shared_ptr<LinearSolver<Matrix, Vector>> coupling_op_solver_;
        std::shared_ptr<Operator<Vector>> S;
        Vector temp;


        bool apply_T(const Vector &rhs, Vector &x)
        {

            assert(!has_nan_or_inf(rhs));

            temp = (*this->B) * rhs;

            assert(!has_nan_or_inf(temp));

            x = local_zeros(local_size(temp));
            bool ok = coupling_op_solver_->apply(temp, x);

            assert(!has_nan_or_inf(x));
            return ok;
        }

        bool apply_T_transpose(const Vector &rhs, Vector &x)
        {
            assert(!has_nan_or_inf(rhs));

            temp = local_zeros(local_size(rhs));
            if(!coupling_op_solver_->apply(rhs, temp)) {
                return false;
            }

            assert(!has_nan_or_inf(temp));

            x = transpose(*this->B) * temp;

            assert(!has_nan_or_inf(x));
            return true;
        }
    };
}


#endif //UTOPIA_SP_STATIC_CONDENSATION_KRYLOV_HPP

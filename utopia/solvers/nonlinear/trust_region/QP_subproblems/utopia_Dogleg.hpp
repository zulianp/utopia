#ifndef UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP
#define UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP
#include "utopia_TRSubproblem.hpp"
#include "utopia_CauchyPoint.hpp"
#include "utopia_Parameters.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

namespace utopia
{

	template<class Matrix, class Vector>
    class Dogleg final: public TRSubproblem<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;
        typedef utopia::LinearSolver<Matrix, Vector>            LinearSolver;

        public:

        Dogleg( const std::shared_ptr <LinearSolver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >()) :
                TRSubproblem<Matrix, Vector>(),
                ls_solver_(linear_solver)
                { }

                inline Dogleg * clone() const override
                {
                    return new Dogleg(std::shared_ptr<LinearSolver>(ls_solver_->clone()));
                }

            bool apply(const Vector &b, Vector &x) override
            {
                return aux_solve(*this->get_operator(), -1.0 * b, x);
            }


        protected:
            bool aux_solve(const Matrix &B, const Vector &g, Vector &p_k) 
            {
                Vector p_N = local_zeros(local_size(p_k)), p_SD = local_zeros(local_size(p_k));
                Scalar g_B_g = dot(g, B * g);

                if(ls_solver_)
                {
                    ls_solver_->solve(B, -1 * g, p_k);
                }
                else
                {
                    utopia_error("Dogleg:: linear solver is missing... \n"); 
                }

                if(norm2(p_k) <= this->current_radius())
                {
                    return true;
                }
                else
                {
                    p_SD = -1.0 *(dot(g, g)/g_B_g) * g;
                    Scalar SD_norm = norm2(p_SD);

                    if(SD_norm > this->current_radius())
                    {
                        Scalar a = SD_norm*SD_norm;
                        Scalar c = - 1.0 * std::pow(this->current_radius(), 2.);

                        Scalar tau = this->quadratic_function(a, 0.0, c);
                        p_k *= tau;

                        return true;
                    }
                    else
                    {
                        Vector d = p_k - p_SD;
                        Scalar d_norm = norm2(d);

                        Scalar a = d_norm*d_norm;
                        Scalar b = 2.0 * dot(p_SD, p_N);
                        Scalar c = SD_norm*SD_norm - 1.0 * std::pow(this->current_radius(), 2.);

                        Scalar tau = this->quadratic_function(a, b, c);
                        p_k = p_SD + tau * d;

                        return true;
                    }
                }
            }

        public:
            void set_linear_solver(const std::shared_ptr<LinearSolver > &ls)
            {
                ls_solver_ = ls;
            }

        private:
            std::shared_ptr<LinearSolver> ls_solver_;
    };
}

#endif //UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP
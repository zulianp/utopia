#ifndef UTOPIA_ROSENBROCK_TRUST_REGION_HPP
#define UTOPIA_ROSENBROCK_TRUST_REGION_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_EigenSolver.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    /**
     * @brief This algorithm is implementation of the following paper:
     * Combining trust region techniques and Rosenbrock methods for gradient systems by Luo, Kelley, Liao, Tam
     */
    template<class Matrix, class Vector>
    class RosenbrockTrustRegion final: public NewtonBase<Matrix, Vector>, public TrustRegionBase<Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef typename NewtonBase<Matrix, Vector>::Solver Solver;
        typedef utopia::EigenSolver<Matrix, Vector>         EigenSolver;

        using NewtonBase<Matrix, Vector>::print_statistics;


    public:
       RosenbrockTrustRegion(   const std::shared_ptr <Solver> &linear_solver,
                                const std::shared_ptr<EigenSolver> & eigen_solver):
                                NewtonBase<Matrix, Vector>(linear_solver),
                                eigen_solver_(eigen_solver),
                                eps_(1e-12),
                                tau_(1e-4)
        {

        }

        Scalar eps() const
        {
            return eps_;
        }

        Scalar tau() const
        {
            return tau_;
        }


        void tau(const Scalar & tau)
        {
            tau_ = tau;
        }

        void eps(const Scalar & eps)
        {
            eps_ = eps;
        }


        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
           using namespace utopia;

            if(this->verbose())
                this->init_solver("RosenbrockTrustRegion", {" it. ", "|| g ||", " E ",  "rho", "lambda", "|| Delta x || "});

            bool converged = false;
            SizeType it = 0;

            Scalar g_norm, s_norm=9e9, H_norm, lambda, lambda_min;
            Scalar ared, pred, rho, expected_reduction;

            Scalar scaling_factor_mat = 1.0 - std::sqrt(2.0)/2.0;
            Scalar scaling_factor_vec = (std::sqrt(2.0)-1.0)/2.0;

            Vector g = local_zeros(local_size(x)), eigenvector_min, s, x_trial;
            Matrix H, H_damped;

            fun.gradient(x, g);
            g_norm = norm2(g);

            Scalar energy_old, energy_new, energy;
            fun.value(x, energy_old);

            fun.hessian(x, H);
            H_norm = norm2(H);

            Matrix I = local_identity(local_size(H).get(0), local_size(H).get(1));

            // tau = 1.0/g_norm;
            // lambda = g_norm;
            Scalar tau = std::min(g_norm, 10.0);
            lambda = 1./tau;

            if(this->verbose())
                PrintInfo::print_iter_status(it, {g_norm, energy_old, 0.0, lambda, 0.0});

            while(!converged)
            {
                // H_damped = (lambda * I);
                // H_damped += (scaling_factor_mat * H);
                H_damped = (scaling_factor_mat * H);
                // H_damped += (lambda * I);
                H_damped.shift_diag(lambda); 

                eigen_solver_->portion_of_spectrum("smallest_real");
                eigen_solver_->number_of_eigenvalues(1);
                eigen_solver_->solve(H_damped);
                eigen_solver_->get_real_eigenpair(0, lambda_min, eigenvector_min);

                if(lambda_min >= eps_)
                {
                    // Rosenbrock step
                    // one could do factorization ones and apply just substitution twice...
                    // but that would assume direct solver but default, which is not really HPC solution...
                    // so for now, we assume cost of one TR iteration is 2 linear solves
                    s = 0 * x;
                    this->linear_solve(H_damped, -1.0 * g, s);
                    Vector x_temp = x + scaling_factor_vec * s;
                    Vector g_temp = local_zeros(local_size(g).get(0));
                    fun.gradient(x_temp, g_temp);
                    s = 0 * x;
                    this->linear_solve(H_damped, -1.0 * g_temp, s);

                    // building trial point
                    x_trial = x + s;

                    pred = -1.0 * dot(g, s) - 0.5 *dot(H* s, s);

                    s_norm = norm2(s);
                    expected_reduction = tau_ * g_norm * std::min(s_norm, g_norm/H_norm);

                    if(pred >= expected_reduction)
                    {
                        fun.value(x_trial, energy_new);
                        ared = energy_old - energy_new;
                        rho = ared/pred;
                    }
                    else
                    {
                        rho = -1.0;
                    }
                }
                else
                {
                    rho = -1.0;
                }

                if(!std::isfinite(rho)){
                    rho = 0.0;
                }

                if(rho > 0.0)
                {
                    x = x_trial;
                    fun.gradient(x, g);
                    g_norm = norm2(g);

                    energy_old = energy_new;
                    energy = energy_new;
                }
                else
                {
                    energy = energy_old;
                }

                // adjusting lambda according to rho
                if(rho < 0.0){
                    lambda *=10.0;
                }
                else if(rho < this->eta1()){
                    lambda *= this->gamma2();
                }
                else if(rho > this->eta2()){
                    lambda *= this->gamma1();
                }

                it++;

                if(this->verbose())
                    PrintInfo::print_iter_status(it, {g_norm, energy,  rho, lambda, s_norm});

                converged = NewtonBase<Matrix, Vector>::check_convergence(it, g_norm, 9e9, s_norm);

                if(!converged && rho > 0.0)
                {
                    fun.hessian(x, H);
                    H_norm = norm2(H);
                }

            } // outer solve loop while(!converged)

            return true;
        }

    void read(Input &in) override
    {
        NewtonBase<Matrix, Vector>::read(in);
        TrustRegionBase<Vector>::read(in);

        in.get("eps", eps_);
        in.get("tau", tau_);

        if(eigen_solver_) {
            in.get("eigen-solver", *eigen_solver_);
        }
    }

    void print_usage(std::ostream &os) const override
    {
        NewtonBase<Matrix, Vector>::print_usage(os);
        TrustRegionBase<Vector>::print_usage(os);

        this->print_param_usage(os, "eps", "real", "Tolerance for checking if smallest eigenvalue.", "1e-12");
        this->print_param_usage(os, "tau", "real", "Constant defining lower bound on predicted reduction.", "1e-12");
    }


    private:
        std::shared_ptr<EigenSolver> eigen_solver_;
        Scalar eps_;
        Scalar tau_; // constant defining lower bound on reduction...


    };

}
#endif //UTOPIA_ROSENBROCK_TRUST_REGION_HPP

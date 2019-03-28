#ifndef UTOPIA_PSEUDO_TRUST_REGION_HPP
#define UTOPIA_PSEUDO_TRUST_REGION_HPP

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
     * Trust region algorithms and timestep selection by D.J. Higham
     */
    template<class Matrix, class Vector>
    class PseudoTrustRegion final: public NewtonBase<Matrix, Vector>, public TrustRegionBase<Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef typename NewtonBase<Matrix, Vector>::Solver Solver;
        typedef utopia::EigenSolver<Matrix, Vector>         EigenSolver;

        using NewtonBase<Matrix, Vector>::print_statistics;


    public:
       PseudoTrustRegion(   const std::shared_ptr <Solver> &linear_solver,
                            const std::shared_ptr<EigenSolver> & eigen_solver):
                            NewtonBase<Matrix, Vector>(linear_solver),
                            eigen_solver_(eigen_solver),
                            eps_(1e-12)
        {

        }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
           using namespace utopia;

            if(this->verbose())
                this->init_solver("PseudoTrustRegion", {" it. ", "|| g ||", " E ",  "rho", "tau", "|| Delta x || "});

            bool converged = false;
            SizeType it = 0;

            Scalar g_norm, s_norm=9e9, tau, lambda_min;
            Scalar ared, pred, rho;

            Vector g = local_zeros(local_size(x)), eigenvector_min, x_trial;
            Vector s = 0 * x;
            Matrix H, H_damped;

            fun.gradient(x, g);
            g_norm = norm2(g);

            Scalar energy_old, energy_new, energy;
            fun.value(x, energy_old);
            fun.hessian(x, H);

            Matrix I = local_identity(local_size(H));

            //tau = 1.0/g_norm;
            tau = std::min(g_norm, 10.0);
            // tau = g_norm;

            if(this->verbose())
                PrintInfo::print_iter_status(it, {g_norm, energy_old, 0.0, tau, 0.0});

            while(!converged)
            {
                H_damped = H + 1./tau * I;

                eigen_solver_->portion_of_spectrum("smallest_real");
                eigen_solver_->number_of_eigenvalues(1);
                eigen_solver_->solve(H_damped);
                eigen_solver_->get_real_eigenpair(0, lambda_min, eigenvector_min);

                if(lambda_min >= eps_)
                {
                    s = 0 * x;
                    this->linear_solve(H_damped, -1.0 * g, s);
                    x_trial = x + s;

                    fun.value(x_trial, energy_new);
                    ared = energy_old - energy_new;

                    pred = -1.0 * dot(g, s) -0.5 *dot(H* s, s);
                    rho = ared/pred;


                    if(rho < this->eta1())
                        tau *= this->gamma1();
                    else if(rho > this->eta2())
                        tau *= this->gamma2();
                }
                else
                {
                    rho = -1.0;
                    tau  *= this->gamma1();
                }

                if(rho > 0.0)
                {
                    x = x_trial;
                    fun.gradient(x, g);
                    energy_old = energy_new;
                    energy = energy_new;

                    norms2(g, s, g_norm, s_norm);
                }
                else
                {
                    s_norm = norm2(s);
                    energy = energy_old;
                }

                it++;

                if(this->verbose())
                {
                    // we can not use s_norm as sometimes is zero, if we do not go for solve...
                    PrintInfo::print_iter_status(it, {g_norm, energy,  rho, tau, s_norm});
                }

                // we can not use s_norm as sometimes is zero, if we do not go for solve...
                converged = NewtonBase<Matrix, Vector>::check_convergence(it, g_norm, 9e9, 9e9);

                if(!converged && rho >0.0){
                    fun.hessian(x, H);
                }

            } // outer solve loop while(!converged)

            return true;
        }


    void read(Input &in) override
    {
        NewtonBase<Matrix, Vector>::read(in);
        in.get("eps", eps_);

        if(eigen_solver_) {
            in.get("eigen-solver", *eigen_solver_);
        }

    }

    void print_usage(std::ostream &os) const override
    {
        NewtonBase<Matrix, Vector>::print_usage(os);
        this->print_param_usage(os, "eps", "real", "Tolerance for checking if smallest eigenvalue.", "1e-12");
    }


    private:
        std::shared_ptr<EigenSolver> eigen_solver_;
        Scalar eps_;


    };

}
#endif //UTOPIA_PSEUDO_TRUST_REGION_HPP

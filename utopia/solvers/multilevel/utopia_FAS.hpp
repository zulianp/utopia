#ifndef UTOPIA_FAS_HPP
#define UTOPIA_FAS_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NewtonBase.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"

#include "utopia_LevelMemory.hpp"


namespace utopia
{
    /**
     * @brief      The class for Full approximation scheme.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class FAS final: public NonlinearMultiLevelBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::NonLinearSmoother<Matrix, Vector>   Smoother;
        typedef utopia::Transfer<Matrix, Vector>   Transfer;

        typedef utopia::NewtonBase<Matrix, Vector>     Solver;
        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;


        public:

        FAS(    const SizeType & n_levels,
                const std::shared_ptr<Smoother> &smoother,
                const std::shared_ptr<Solver> &coarse_solver):
                NonlinearMultiLevelBase<Matrix,Vector>(n_levels),
                smoother_(smoother),
                coarse_solver_(coarse_solver)
        {

        }


        std::string name() override
        {
            return "FAS";
        }

        void read(Input &in) override
        {
            NonlinearMultiLevelBase<Matrix, Vector>::read(in);

            if(smoother_)
            {
                in.get("smoother", *smoother_);
            }
            if(coarse_solver_)
            {
                in.get("coarse_solver", *coarse_solver_);
            }
        }

        void print_usage(std::ostream &os) const override
        {
            NonlinearMultiLevelBase<Matrix, Vector>::print_usage(os);

            this->print_param_usage(os, "coarse_solver", "NewtonBase", "Input parameters for coarse level QP solvers.", "-");
            this->print_param_usage(os, "smoother", "NonLinearSmoother", "Input parameters for fine level QP solver.", "-");
        }


        bool solve(Vector & x_h) override
        {
            bool converged = false;
            SizeType it = 0, n_levels = this->n_levels();
            Scalar r_norm, r0_norm=1, rel_norm=1, energy;

            if(this->transfers_.size() + 1 != this->level_functions_.size()){
                utopia_error("FAS::solve size of transfer and level functions do not match... \n");
            }

            std::string header_message = this->name() + ": " + std::to_string(n_levels) +  " levels";
            this->init_solver(header_message, {" it. ", "|| grad ||", "r_norm" , "Energy"});

            SizeType loc_size = local_size(x_h).get(0);
            this->init_memory(loc_size);


            memory_.x[n_levels-1] = x_h;
            memory_.g[n_levels-1] = local_zeros(local_size(memory_.x[n_levels-1]));


            this->function(n_levels-1).gradient(memory_.x[n_levels-1], memory_.g[n_levels-1]);
            r0_norm = norm2(memory_.g[n_levels-1]);
            r_norm = r0_norm;

            this->function(n_levels-1).value(memory_.x[n_levels-1], energy);

            if(this->verbose())
                PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy});

            it++;

            while(!converged)
            {
                this->multiplicative_cycle(n_levels);

#ifdef CHECK_NUM_PRECISION_mode
                if(has_nan_or_inf(memory_.x[n_levels-1]) == 1)
                {
                    memory_.x[n_levels-1] = local_zeros(local_size(memory_.x[n_levels-1]));
                    return true;
                }
#endif

                this->function(n_levels-1).gradient(memory_.x[n_levels-1], memory_.g[n_levels-1]);
                this->function(n_levels-1).value(memory_.x[n_levels-1], energy);

                r_norm = norm2(memory_.g[n_levels-1]);
                rel_norm = r_norm/r0_norm;

                // print iteration status on every iteration
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy});

                // check convergence and print interation info
                converged = this->check_convergence(it, r_norm, rel_norm, 1);
                it++;
            }

            this->print_statistics(it);

#ifdef CHECK_NUM_PRECISION_mode
            if(has_nan_or_inf(memory_.x[n_levels-1]) == 1)
                exit(0);
#endif
            x_h = memory_.x[n_levels-1];

            return true;
        }



    protected:

        void init_memory(const SizeType & fine_local_size) override
        {
            memory_.init(this->n_levels());
            memory_.g_diff[this->n_levels()-1] = local_zeros(fine_local_size);
        }


        bool multiplicative_cycle(const SizeType & l)
        {
            for(auto l = this->n_levels()-1; l > 0; l--)
            {
                // pre-smoothing
                smoothing(this->function(l), memory_.x[l], memory_.g_diff[l], this->pre_smoothing_steps());

                // building intial guess
                this->transfer(l-1).project_down(memory_.x[l], memory_.x[l-1]);
                memory_.x_0[l-1] = memory_.x[l-1];

                // multilevel gradient ...
                this->function(l).gradient(memory_.x[l], memory_.g[l]);

                memory_.g[l] -= memory_.g_diff[l];

                this->transfer(l-1).restrict(memory_.g[l], memory_.g_diff[l-1]);
                this->function(l-1).gradient(memory_.x[l-1], memory_.g[l-1]);

                memory_.g_diff[l-1] = memory_.g[l-1] - memory_.g_diff[l-1];

                this->zero_correction_related_to_equality_constrain(this->function(l-1), memory_.g_diff[l-1]);
            }

            coarse_solve(this->function(0), memory_.x[0],memory_.g_diff[0]);

            for(auto l = 0; l < this->n_levels()-1; l++)
            {
                memory_.c[l] = memory_.x[l] - memory_.x_0[l];
                this->transfer(l).interpolate(memory_.c[l], memory_.c[l+1]);

                this->zero_correction_related_to_equality_constrain(this->function(l+1), memory_.c[l+1]);

                memory_.x[l+1] += memory_.c[l+1];
                smoothing(this->function(l+1), memory_.x[l+1], memory_.g_diff[l+1], this->pre_smoothing_steps());
            }

            return true;
        }

        bool smoothing(Function<Matrix, Vector> &fun,  Vector &x, const Vector &f, const SizeType & nu = 1)
        {
            smoother_->sweeps(nu);
            smoother_->smooth(fun, x, f);
            return true;
        }

        bool coarse_solve(Fun &fun, Vector &x, const Vector & rhs)
        {
            coarse_solver_->max_it(1);
            coarse_solver_->solve(fun, x, rhs);
            return true;
        }

    public:
        bool change_coarse_solver(const std::shared_ptr<Solver> &nonlinear_solver = std::shared_ptr<Solver>())
        {
            coarse_solver_ = nonlinear_solver;
            return true;
        }



    protected:
        std::shared_ptr<Smoother>           smoother_;
        std::shared_ptr<Solver>             coarse_solver_;

    private:
        FASLevelMemory <Matrix, Vector>         memory_;


    };

}

#endif //UTOPIA_FAS_HPP


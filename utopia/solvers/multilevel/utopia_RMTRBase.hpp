#ifndef UTOPIA_RMTR_BASE_HPP
#define UTOPIA_RMTR_BASE_HPP

#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"

#include "utopia_TRSubproblem.hpp"
#include "utopia_Linear.hpp"
#include "utopia_Level.hpp"
#include "utopia_LS_Strategy.hpp"

#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_TRBase.hpp"

#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_LevelMemory.hpp"

namespace utopia
{
    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class RMTRBase : public NonlinearMultiLevelBase<Matrix, Vector>
    {
        public:
            typedef UTOPIA_SCALAR(Vector)                       Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;
            typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

        RMTRBase(const SizeType & n_levels): NonlinearMultiLevelBase<Matrix,Vector>(n_levels),
                                            _it_global(0),
                                            _skip_BC_checks(false),
                                            _max_coarse_it(2),
                                            _max_sucessful_smoothing_it(1),
                                            _max_sucessful_coarse_it(1),
                                            _max_QP_smoothing_it(5),
                                            _max_QP_coarse_it(50),
                                            _eps_delta_termination(0.001),
                                            _grad_smoothess_termination(0.5),
                                            _check_gradient_smoothness(true),
                                            _hessian_update_delta(0.15),
                                            _hessian_update_eta(0.5),
                                            _verbosity_level(VERBOSITY_LEVEL_NORMAL),
                                            _norm_schedule(ALWAYS), 
                                            red_(FG_LIGHT_MAGENTA),
                                            def_(FG_DEFAULT),
                                            yellow_(FG_LIGHT_YELLOW),
                                            green_(FG_LIGHT_GREEN), 
                                            ml_derivs_(n_levels)
        {

        }


        virtual ~RMTRBase(){}

        virtual void read(Input &in) override
        {
            NonlinearMultiLevelBase<Matrix, Vector>::read(in);

            in.get("skip_BC_checks", _skip_BC_checks);
            in.get("max_coarse_it", _max_coarse_it);
            in.get("max_sucessful_smoothing_it", _max_sucessful_smoothing_it);
            in.get("max_sucessful_coarse_it", _max_sucessful_coarse_it);
            in.get("max_QP_smoothing_it", _max_QP_smoothing_it);
            in.get("max_QP_coarse_it", _max_QP_coarse_it);

            in.get("eps_delta_termination", _eps_delta_termination);
            in.get("grad_smoothess_termination", _grad_smoothess_termination);
            in.get("check_gradient_smoothness", _check_gradient_smoothness);

            in.get("hessian_update_delta", _hessian_update_delta);
            in.get("hessian_update_eta", _hessian_update_eta);
        }

        virtual void print_usage(std::ostream &os) const override
        {
            NonlinearMultiLevelBase<Matrix, Vector>::print_usage(os);

            this->print_param_usage(os, "skip_BC_checks", "bool", "Skip treatment of BC conditions.", "false");
            this->print_param_usage(os, "max_coarse_it", "int", "Maximum number of coarse iterations.", "2");
            this->print_param_usage(os, "max_sucessful_smoothing_it", "int", "Maximum number of smoothing steps.", "1");
            this->print_param_usage(os, "max_sucessful_coarse_it", "int", "Maximum number of corse iterations.", "1");

            this->print_param_usage(os, "max_QP_smoothing_it", "int", "Maximum number of iterations for fine level QP solver.", "5");
            this->print_param_usage(os, "max_QP_coarse_it", "int", "Maximum number of iterations for coarse level QP solver.", "50");

            this->print_param_usage(os, "eps_delta_termination", "real", "Constant used for quiting recursion based on size of tr. radius.", "0.001");
            this->print_param_usage(os, "grad_smoothess_termination", "real", "Constant used for quiting recursion based on smoothness of the gradient on given level.", "0.5");
            this->print_param_usage(os, "check_gradient_smoothness", "real", "Flag turning on/off check for gradient smoothiness.", "true");

            this->print_param_usage(os, "hessian_update_delta", "real", "Constant used for deciding whether to assemble fresh hessian or no.", "0.15");
            this->print_param_usage(os, "hessian_update_eta", "real", "Constant used for deciding whether to assemble fresh hessian or no.", "0.5");
            this->print_param_usage(os, "verbosity_level", "VerbosityLevel", "Specifies level of verbosity.", "VERBOSITY_LEVEL_NORMAL");
        }


        virtual VerbosityLevel verbosity_level() const
        {
            return _verbosity_level;
        }

        virtual void verbosity_level(const VerbosityLevel & level )
        {
            _verbosity_level = this->verbose() ? level : VERBOSITY_LEVEL_QUIET;
        }

        virtual MultilevelNormSchedule norm_schedule() const
        {
            return _norm_schedule;
        }

        virtual void norm_schedule(const MultilevelNormSchedule & schedule)
        {
            _norm_schedule = schedule; 
        }        

        virtual void set_grad_smoothess_termination(const Scalar & grad_smoothess_termination)
        {
            _grad_smoothess_termination = grad_smoothess_termination;
        }

        virtual Scalar  get_grad_smoothess_termination( ) const
        {
            return _grad_smoothess_termination;
        }


        virtual void skip_BC_checks(const bool skip_BC_checks)
        {
            _skip_BC_checks = skip_BC_checks;
        }

        virtual bool skip_BC_checks() const
        {
            return _skip_BC_checks;
        }

        virtual void max_coarse_it(const SizeType & max_coarse_it)
        {
            _max_coarse_it = max_coarse_it;
        }

        virtual SizeType max_coarse_it() const
        {
            return _max_coarse_it;
        }



        // this parameter define total number of smoothing its - sucessful and unsuccessful
        // while pre_smooting and post_smoothing count just sucessful iterations
        virtual void max_sucessful_smoothing_it(const SizeType & max_sucessful_smoothing_it)
        {
            _max_sucessful_smoothing_it = max_sucessful_smoothing_it;
        }


        virtual SizeType max_sucessful_smoothing_it() const
        {
            return _max_sucessful_smoothing_it;
        }

        virtual void max_sucessful_coarse_it(const SizeType & max_sucessful_coarse_it)
        {
            _max_sucessful_coarse_it = max_sucessful_coarse_it;
        }        

        virtual SizeType max_sucessful_coarse_it() const
        {
            return _max_sucessful_coarse_it;
        }        


        virtual void max_QP_smoothing_it(const SizeType & num_it)
        {
            _max_QP_smoothing_it = num_it;
        }


        virtual SizeType max_QP_smoothing_it() const
        {
            return _max_QP_smoothing_it;
        }        

        virtual void max_QP_coarse_it(const SizeType & num_it)
        {
            _max_QP_coarse_it = num_it;
        }

        virtual SizeType max_QP_coarse_it() const
        {
            return _max_QP_coarse_it;
        }


        virtual void check_grad_smoothness(const bool flg)
        {
            _check_gradient_smoothness = flg;
        }

        virtual bool check_grad_smoothness() const
        {
            return _check_gradient_smoothness;
        }


        virtual void eps_delta_termination(const Scalar &eps)
        {
            _eps_delta_termination = eps;
        }

        virtual Scalar eps_delta_termination() const
        {
            return _eps_delta_termination;
        }

        virtual void hessian_update_delta(const Scalar &delta)
        {
            _hessian_update_delta = delta;
        }

        virtual Scalar hessian_update_delta() const
        {
            return _hessian_update_delta;
        }

        virtual void hessian_update_eta(const Scalar &eps)
        {
            _hessian_update_eta = eps;
        }

        virtual Scalar hessian_update_eta() const
        {
            return _hessian_update_eta;
        }                



    protected:
        virtual bool get_multilevel_hessian(const Fun & fun, const SizeType & level)
        {
            return  ml_derivs_.compute_hessian(level, fun, memory_.x[level]);
        }


        virtual bool get_multilevel_gradient(const Fun & fun, const Vector & s_global, const SizeType & level)
        {
            return ml_derivs_.compute_gradient(level, fun, memory_.x[level], s_global);
        }


        virtual Scalar get_multilevel_energy(const Fun & fun, const Vector & s_global, const SizeType & level)
        {
            return ml_derivs_.compute_energy(level, fun, memory_.x[level], s_global);
        }


        virtual void init_memory(const SizeType & /*fine_local_size */)
        {
            const std::vector<SizeType> & dofs =  this->local_level_dofs(); 

            ml_derivs_.init_memory(dofs); 
            memory_.init_memory(dofs); 
        }



        virtual void print_level_info(const SizeType & level)
        {
            if(verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
            {
                if(level == 0)
                {
                    std::cout << yellow_;
                    std::string solver_type = "COARSE SOLVE:: " + std::to_string(level);
                    this->print_init_message(solver_type, {" it. ", "|| g ||", "   E + <g_diff, x>", "ared   ",  "  pred  ", "  rho  ", "  delta "});
                }
                else
                {
                    std::cout << green_;
                    std::string solver_type = "SMOOTHER:  " + std::to_string(level);
                    this->print_init_message(solver_type, {" it. ", "|| g ||", "   E + <g_diff, x>", "ared   ",  "  pred  ", "  rho  ", "  delta "});
                }
            }
        }


    protected:
        SizeType                         _it_global;                 /** * global iterate counter  */
        

        bool                            _skip_BC_checks;

        SizeType                        _max_coarse_it;             /** * maximum iterations on coarse level   */
        SizeType                        _max_sucessful_smoothing_it;          /** * max smoothing iterations  */
        SizeType                        _max_sucessful_coarse_it;          /** * max smoothing iterations  */

        SizeType                        _max_QP_smoothing_it;
        SizeType                        _max_QP_coarse_it;


        Scalar                         _eps_delta_termination;      /** * maximum delta allowed on coarse level - makes sure that coarse level corection stays inside fine level radius  */

        Scalar                         _grad_smoothess_termination; /** * determines when gradient is not smooth enough => does pay off to go to coarse level at all  */
        bool                           _check_gradient_smoothness;

        Scalar                         _hessian_update_delta;       /** * tolerance used for updating hessians */
        Scalar                         _hessian_update_eta;         /** * tolerance used for updating hessians */

        VerbosityLevel                  _verbosity_level;
        MultilevelNormSchedule          _norm_schedule; 

        ColorModifier red_;
        ColorModifier def_;
        ColorModifier yellow_;
        ColorModifier green_;

        RMTRLevelMemory<Matrix, Vector>         memory_;
        MultilevelDerivEval<Matrix, Vector, CONSISTENCY_LEVEL>  ml_derivs_; 

    };

}

#endif //UTOPIA_RMTR_BASE_HPP
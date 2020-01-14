// #ifndef UTOPIA_QUASI_RMTR_INF_HPP
// #define UTOPIA_QUASI_RMTR_INF_HPP

// #include "utopia_NonLinearSmoother.hpp"
// #include "utopia_NonLinearSolver.hpp"
// #include "utopia_Core.hpp"
// #include "utopia_NonlinearMultiLevelBase.hpp"

// #include "utopia_TRSubproblem.hpp"
// #include "utopia_TrustRegionVariableBound.hpp"

// #include "utopia_Linear.hpp"
// #include "utopia_Level.hpp"

// #include "utopia_NonLinearSolver.hpp"
// #include "utopia_NonLinearSmoother.hpp"
// #include "utopia_TRBase.hpp"

// #include "utopia_MultiLevelEvaluations.hpp"
// #include "utopia_RMTR.hpp"

// #include "utopia_Transfer.hpp"
// #include "utopia_IdentityTransfer.hpp"


// namespace utopia
// {
//     /**
//      * @brief      The class for Quasi RMTR in infinity norm...
//      *
//      * @tparam     Matrix
//      * @tparam     Vector
//      */
//     template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
//     class QuasiRMTR_inf :   public RMTR<Matrix, Vector, CONSISTENCY_LEVEL>,  public MultilevelVariableBoundSolverInterface<Matrix, Vector>
//     {
//         typedef UTOPIA_SCALAR(Vector)                       Scalar;
//         typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

//         typedef utopia::MatrixFreeQPSolver<Vector>          TRSubproblem;
//         typedef std::shared_ptr<TRSubproblem>               TRSubproblemPtr;

//         typedef utopia::Transfer<Matrix, Vector>            Transfer;
//         typedef utopia::Level<Matrix, Vector>               Level;

//         typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

//         typedef utopia::BoxConstraints<Vector>                          BoxConstraints;
//         typedef utopia::RMTR<Matrix, Vector, CONSISTENCY_LEVEL>         RMTR;

//         typedef utopia::HessianApproximation<Vector>    HessianApproximation;
//         typedef std::shared_ptr<HessianApproximation>   HessianApproxPtr;

//         typedef utopia::MultilevelVariableBoundSolverInterface<Matrix, Vector>  MLConstraints;

//         static_assert(utopia::is_first_order<CONSISTENCY_LEVEL>::value, "utopia::QuasiRMTR does not support second order, nor galerkin consistency, nor Galerkin.");

//     public:


//         QuasiRMTR_inf(  const SizeType & n_levels):
//                     RMTR(n_levels)
//         {
//             hessian_approxs_.resize(n_levels);
//         }

//         virtual ~QuasiRMTR_inf()
//         {
//             // do we need to destroy some memory or no???
//         }


//         virtual void read(Input &in) override
//         {
//             RMTR::read(in);

//             if(_tr_subproblems.size() > 0)
//             {
//                 in.get("coarse-QPSolver", *_tr_subproblems[0]);

//                 for(auto i=1; i < _tr_subproblems.size(); i++)
//                     in.get("fine-QPSolver", *_tr_subproblems[i]);
//             }

//             if(hessian_approxs_.size() > 0)
//             {
//                 for(auto i=1; i < hessian_approxs_.size(); i++)
//                     in.get("hessian-approx-strategy", *hessian_approxs_[i]);
//             }
//         }

//         virtual void print_usage(std::ostream &os) const override
//         {
//             RMTR::print_usage(os);

//             this->print_param_usage(os, "coarse-QPSolver", "MatrixFreeQPSolver", "Input parameters for fine level QP solvers.", "-");
//             this->print_param_usage(os, "fine-QPSolver", "MatrixFreeQPSolver", "Input parameters for coarse level QP solver.", "-");
//             this->print_param_usage(os, "hessian-approx-strategy", "HessianApproximation", "Input parameters for hessian approximation strategies.", "-");
//         }


//         virtual std::string name() override
//         {
//             return "QuasiRMTR_inf";
//         }

//         virtual bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy)
//         {
//             if(hessian_approxs_.size() != this->n_levels())
//                 hessian_approxs_.resize(this->n_levels());

//             for(auto l  = 0; l != hessian_approxs_.size(); ++l)
//                 hessian_approxs_[l] = std::shared_ptr<HessianApproximation>(strategy->clone());

//             return true;
//         }

//         virtual bool set_hessian_approximation_strategies(const std::vector<HessianApproxPtr> &strategies)
//         {
//             if(strategies.size() != this->n_levels()){
//                 utopia_error("utopia::QuasiRMTR::set_hessian_approximation_strategies:: Number of strategies does not equal with levels in ML hierarchy. \n");
//             }

//             hessian_approxs_ = strategies;

//             return true;
//         }

//         virtual bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy, const SizeType & level)
//         {
//             if(hessian_approxs_.size() != this->n_levels())
//                 hessian_approxs_.resize(this->n_levels());

//             if(level <= this->n_levels())
//             {
//                 hessian_approxs_[level] = strategy;
//             }
//             else{
//                 utopia_error("utopia::QuasiRMTR::set_tr_strategy:: Requested level exceeds number of levels in ML hierarchy. \n");
//             }

//             return true;
//         }



//         bool set_coarse_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy)
//         {
//             if(_tr_subproblems.size() != this->n_levels())
//                 _tr_subproblems.resize(this->n_levels());

//             _tr_subproblems[0] = strategy;

//             return true;
//         }

//         bool set_fine_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy)
//         {
//             if(_tr_subproblems.size() != this->n_levels())
//                 _tr_subproblems.resize(this->n_levels());

//             // starting from level 1 ....
//             for(std::size_t l = 1; l != _tr_subproblems.size(); ++l)
//                 _tr_subproblems[l] = std::shared_ptr<TRSubproblem>(strategy->clone());


//             return true;
//         }


//         bool set_tr_strategies(const std::vector<TRSubproblemPtr> &strategies)
//         {
//             if(strategies.size() != this->n_levels()){
//                 utopia_error("utopia::RMTR::set_tr_strategies:: Number of tr strategies MUST be equal to number of levels in ML hierarchy. \n");
//             }

//             _tr_subproblems = strategies;

//             return true;
//         }



//     protected:
//         virtual void init_memory() override
//         {
//             RMTR::init_memory();

//             // TODO:: modify allocation of constraints 
//             const std::vector<SizeType> & dofs =  this->local_level_dofs(); 
//             MLConstraints::init_constr_memory(this->n_levels(), dofs.back()); 

//             const SizeType fine_level = this->n_levels()-1;

//             // precompute norms of prolongation operators needed for projections of constraints...
//             for(auto l = 0; l < fine_level; l++){
//                 this->constraints_memory_.P_inf_norm[l] = this->transfer(l).interpolation_inf_norm();
//             }

//             // Hessian approximations 
//             std::cout<< "RMTR init memory ---"<< fine_level << " ---  "<< this->memory_.x[fine_level].size() << "  \n"; 

//             hessian_approxs_[fine_level]->initialize(this->memory_.x[fine_level], this->memory_.g[fine_level]);
//         }



//         // since TR bounds are weak bounds...
//         virtual bool check_feasibility(const SizeType & level ) override
//         {
//             return MLConstraints::check_feasibility(level, this->memory_.x[level]); 
//         }


//         virtual void init_level(const SizeType & level) override
//         {
//             RMTR::init_level(level);

//             const SizeType finer_level = level+1;
//             MLConstraints::init_level(level, this->memory_.x[finer_level], this->memory_.x[level], this->memory_.delta[finer_level]); 
//         }


//         virtual bool get_multilevel_hessian(const Fun & /*fun*/, const SizeType & /*level*/) override
//         {
//             return false;
//         }

//         // -------------------------- tr radius managment ---------------------------------------------
//         /**
//          * @brief      Updates delta on given level
//          *
//          * @param[in]  rho        The rho
//          * @param[in]  level      The level
//          * @param[in]  s_global   Sum of all corrections on given level
//          * @param      converged  convergence flag
//          */
//         virtual bool delta_update(const Scalar & rho, const SizeType & level, const Vector & /*s_global*/) override
//         {
//             Scalar intermediate_delta;

//             // we could do also more sophisticated options, but lets not care for the moment ...
//             if(rho < this->eta1())
//                  intermediate_delta = std::max(this->gamma1() * this->memory_.delta[level], 1e-15);
//             else if (rho > this->eta2() )
//                  intermediate_delta = std::min(this->gamma2() * this->memory_.delta[level], 1e15);
//             else
//                 intermediate_delta = this->memory_.delta[level];

//             this->memory_.delta[level] = intermediate_delta;

//             return false;
//         }


//         /**
//          * @brief      "Heuristics", which decides if it makes sense to go to the coarse level or no
//          *
//          * @param[in]  g_restricted  Restricted gradient
//          * @param[in]  g_coarse      Coarse level gradient
//          *
//          */
//         virtual bool recursion_termination_smoothness(const Vector & g_restricted, const Vector & g_coarse, const SizeType & level) override
//         {
//             //FIXME remove temporary
//             Vector Pc;

//             Vector x_g = this->memory_.x[level] - g_restricted;
//             MLConstraints::get_projection(x_g, this->constraints_memory_.active_lower[level], this->constraints_memory_.active_upper[level], Pc);
//             Pc -= this->memory_.x[level];
//             Scalar Rg_norm =  norm2(Pc);


//             x_g = this->memory_.x[level] - g_coarse;
//             MLConstraints::get_projection(x_g, this->constraints_memory_.active_lower[level], this->constraints_memory_.active_upper[level], Pc);
//             Pc -= this->memory_.x[level];
//             Scalar  g_norm =  norm2(Pc);

//             return (Rg_norm >= this->grad_smoothess_termination() * g_norm) ? true : false;
//         }


//         // measuring wrt to feasible set...
//         virtual Scalar criticality_measure(const SizeType & level) override
//         {
//             return MLConstraints::criticality_measure_inf(level, this->memory_.x[level], this->memory_.g[level]); 
//         }

//         // TODO:: verify + if correct grad_smoothess_termination from above can be removed 
//         virtual void criticality_measures(const SizeType & level, const Vector & g_restricted, const Vector & g_coarse, Scalar & Rg_norm, Scalar & g_norm) override
//         {
//             // norms2(g_restricted, g_coarse, Rg_norm, g_norm);

//             Rg_norm =  MLConstraints::criticality_measure_inf(level, this->memory_.x[level], g_restricted); 
//             g_norm  =  MLConstraints::criticality_measure_inf(level, this->memory_.x[level], g_coarse); 
//         }         


//         /**
//          * @brief      Solves TR subroblem for given level
//          *
//          * @param[in]  H      The hessian
//          * @param[in]  g      The gradient
//          * @param      s      New correction
//          * @param[in]  level  The level
//          *
//          */
//         virtual bool solve_qp_subproblem(const SizeType & level, const bool & flg) override
//         {
//             Scalar radius = this->memory_.delta[level];

//             // first we need to prepare box of intersection of level constraints with tr. constraints
//             Vector l = this->constraints_memory_.active_lower[level] - this->memory_.x[level];
//             each_transform(l, l, [radius](const SizeType /*i*/, const Scalar val) -> Scalar {
//                 return (val >= -1*radius)  ? val : -1 * radius;  }
//             );

//             Vector u =  this->constraints_memory_.active_upper[level] - this->memory_.x[level];
//             each_transform(u, u, [radius](const SizeType /*i*/, const Scalar val) -> Scalar {
//               return (val <= radius)  ? val : radius; }
//             );


//             //FIXME remove the fact that you are creating copies of l and u
//             // generating constraints to go for QP solve
//             auto box = make_box_constaints(std::make_shared<Vector>(l), std::make_shared<Vector>(u));


//             // cast to iterative solver to do this... 
//             // this->_tr_subproblems[level]->atol(1e-14);

//             // //To do this, we need to do some casting to QP solver, not to matrix free thing... 
//             // if(flg){
//             //     this->_tr_subproblems[level]->max_it(this->max_QP_coarse_it());
//             // }
//             // else{
//             //     this->_tr_subproblems[level]->max_it(this->max_QP_smoothing_it());
//             // }


//             auto multiplication_action = hessian_approxs_[level]->build_apply_H();
//             _tr_subproblems[level]->set_box_constraints(box);

//             this->_tr_subproblems[level]->solve(*multiplication_action, -1.0 * this->memory_.g[level], this->memory_.s[level]);


//             // ----- just for debugging pourposes ----------------
//             Vector s_old = this->memory_.s[level]; 
//             MLConstraints::get_projection(s_old, l, u, this->memory_.s[level]); 



//             return true;
//         }


//         virtual Scalar get_pred(const SizeType & level) override
//         {
//             Scalar l_term = dot(this->memory_.g[level], this->memory_.s[level]);
//             Scalar qp_term = hessian_approxs_[level]->compute_uHu_dot(this->memory_.s[level]);
//             return  (- l_term - 0.5 * qp_term);
//         }


//         virtual bool update_level(const SizeType & level) override
//         {
//             Vector grad_old = this->memory_.g[level];
//             this->get_multilevel_gradient(this->function(level), this->memory_.s_working[level], level);
//             Vector y = this->memory_.g[level] - grad_old;

//             // swap back....
//             this->memory_.g[level] = grad_old;

//             hessian_approxs_[level]->update(this->memory_.s[level], y, this->memory_.x[level], this->memory_.g[level]);

//             return true;
//         }


//         virtual void initialize_local_solve(const SizeType & level, const LocalSolveType & solve_type) override
//         {

//             // std::cout<<"initialize_local_solve  \n"; 

//             // if(!(solve_type == PRE_SMOOTHING && level == this->n_levels()-1))
//             // {
//             //     // this is interesting heuristic
//             //     if(solve_type == PRE_SMOOTHING || solve_type == COARSE_SOLVE)
//             //     {
//                     hessian_approxs_[level]->reset();
//                     // hessian_approxs_[level]->initialize();

//                     // std::cout<<"init local solve: "<< level << " sizee:  "<< this->memory_.x[level].size(); 

//                     hessian_approxs_[level]->initialize(this->memory_.x[level], this->memory_.g[level]);    
//             //     }
//             // }
//         }

//         virtual bool check_initialization() override
//         {
//             RMTR::check_initialization();

//             if(this->hessian_approxs_.size() != this->n_levels()){
//                 utopia_error("utopia::QuasiRMTR:: number of Hessian approximation streategies and levels not equal. \n");
//                 return false;
//             }

//             return true;
//         }


//     private:
//         std::vector<TRSubproblemPtr>            _tr_subproblems;

//     protected:
//         std::vector<HessianApproxPtr>  hessian_approxs_;

//     };

// }

// #endif // UTOPIA_QUASI_RMTR_INF_HPP
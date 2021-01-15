#ifndef UTOPIA_RMTR_BASE_IMPL_HPP
#define UTOPIA_RMTR_BASE_IMPL_HPP

#include "utopia_Backtracking.hpp"
#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_Function.hpp"
#include "utopia_LevelMemory.hpp"

namespace utopia {
    template <class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL>
    bool RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>::solve(Vector &x_h) {
        if (!this->check_initialization()) {
            return false;
        }

        bool converged = false;
        SizeType fine_level = this->n_levels() - 1;

        //-------------- INITIALIZATIONS ---------------
        this->init_memory();
        this->memory_.x[fine_level] = x_h;

        if (!this->skip_BC_checks()) {
            this->make_iterate_feasible(this->function(fine_level), this->memory_.x[fine_level]);
            this->handle_equality_constraints();
        }

        this->memory_.energy[fine_level] =
            this->get_multilevel_gradient_energy(this->function(fine_level), fine_level, this->memory_.s[fine_level]);
        this->memory_.gnorm[fine_level] = this->criticality_measure(fine_level);
        this->_it_global = 0;

        //----------------------------------------------
        if (this->verbosity_level() >= VERBOSITY_LEVEL_NORMAL && this->verbose() == true && mpi_world_rank() == 0) {
            std::cout << this->red_;
            std::string name_id =
                this->name() + "     Number of levels: " + std::to_string(fine_level + 1) +
                "   \n Fine level local dofs: " + std::to_string(this->local_level_layouts_.back().local_size()) + "/" +
                std::to_string(x_h.size());
            this->init_solver(name_id, {" it. ", "|| g ||", "   E "});

            PrintInfo::print_iter_status(this->_it_global,
                                         {this->memory_.gnorm[fine_level], this->memory_.energy[fine_level]});
            std::cout << this->def_;
        }

        UTOPIA_NO_ALLOC_BEGIN("RMTR::region1");
        while (!converged) {
            if (this->cycle_type() == MULTIPLICATIVE_CYCLE)
                this->multiplicative_cycle(fine_level);
            else {
                std::cout << "ERROR::UTOPIA_RMTR << unknown cycle type, solving in "
                             "multiplicative manner ... \n";
                this->multiplicative_cycle(fine_level);
            }

#ifdef CHECK_NUM_PRECISION_mode
            if (has_nan_or_inf(this->memory_.x[fine_level]) == true) {
                this->memory_.x[fine_level].set(0.0);
                return false;
            }
#endif

            this->_it_global++;

            if (this->verbose() && mpi_world_rank() == 0) {
                std::cout << this->red_;
                if (this->verbosity_level() > VERBOSITY_LEVEL_NORMAL && this->verbose() == true) {
                    this->print_init_message("RMTR OUTER SOLVE", {" it. ", "|| g ||", "   E "});
                }

                PrintInfo::print_iter_status(this->_it_global,
                                             {this->memory_.gnorm[fine_level], this->memory_.energy[fine_level]});
                std::cout << this->def_;
            }

            // check convergence
            converged = this->check_global_convergence(
                this->_it_global, this->memory_.gnorm[fine_level], 9e9, this->memory_.delta[fine_level]);
        }

        UTOPIA_NO_ALLOC_END();

        // benchmarking
        NonlinearMultiLevelBase<Matrix, Vector>::print_statistics(this->_it_global);

        // passing result outside
        x_h = this->memory_.x[fine_level];
        return true;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template <class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL>
    bool RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>::multiplicative_cycle(const SizeType &level) {
        Scalar ared = 0.0, coarse_reduction = 0.0, rho = 0.0;
        Scalar E_old = 0.0, E_new = 0.0;
        bool converged = false, smoothness_flg = true;

        //----------------------------------------------------------------------------
        //                   presmoothing
        //----------------------------------------------------------------------------

        if (this->pre_smoothing_steps() != 0) {
            UTOPIA_NO_ALLOC_BEGIN("RMTR::region2");
            converged = this->local_tr_solve(level, PRE_SMOOTHING);
            UTOPIA_NO_ALLOC_END();
        } else {
            converged = false;
        }

        // making sure that correction does not exceed tr radius ...
        if (converged) {
            if (this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE) {
                std::cout << "level converged after pre-smoothing \n";
            }

            return true;
        }

        if (this->pre_smoothing_steps() == 0 && level < this->n_levels() - 1) {
            // s_global is assumed to be zero at this point
            this->get_multilevel_gradient(this->function(level), level);
            this->memory_.gnorm[level] = this->criticality_measure(level);
        }

        if (level == this->n_levels() - 1) {
            converged = this->criticality_measure_termination(level);
            if (converged == true) {
                if (this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE) {
                    std::cout << "level converged after pre-smoothing, "
                                 "criticality_measure_termination \n";
                }
                return true;
            }
        }

        UTOPIA_NO_ALLOC_BEGIN("RMTR::region11");

        smoothness_flg = this->init_consistency_terms(level, this->memory_.energy[level]);

        UTOPIA_NO_ALLOC_END();

        // UTOPIA_NO_ALLOC_BEGIN("RMTR::region12");
        //----------------------------------------------------------------------------
        //                   additional coarse level initialization...
        //----------------------------------------------------------------------------
        this->memory_.x_0[level - 1] = this->memory_.x[level - 1];

        // at this point s_global on coarse level is 0, so we can simplify - NOT TRUE
        // IF MG_OPT TYPE OF ML_MODEL IS USED -> better to do directly in specialized
        // class
        // coarse_reduction =
        //     this->get_multilevel_energy(this->function(level - 1), level - 1);

        coarse_reduction = this->memory_.energy[level - 1];

        // store energy in order to avoid evaluation in the first local_solve
        this->memory_.energy[level - 1] = coarse_reduction;
        // UTOPIA_NO_ALLOC_END();

        //----------------------------------------------------------------------------
        //               recursion  / Taylor correction
        //----------------------------------------------------------------------------

        if (level == 1 && smoothness_flg) {
            this->local_tr_solve(level - 1, COARSE_SOLVE);
        } else if (smoothness_flg) {
            // recursive call into RMTR
            for (SizeType k = 0; k < this->mg_type(); k++) {
                SizeType l_new = level - 1;
                this->multiplicative_cycle(l_new);
            }
        }

        if (smoothness_flg) {
            //----------------------------------------------------------------------------
            //                       building trial point from coarse level
            //----------------------------------------------------------------------------
            // UTOPIA_NO_ALLOC_BEGIN("RMTR::region13");
            this->memory_.s[level - 1] = this->memory_.x[level - 1] - this->memory_.x_0[level - 1];
            coarse_reduction -= this->memory_.energy[level - 1];

            this->transfer(level - 1).interpolate(this->memory_.s[level - 1], this->memory_.s[level]);

            if (!this->skip_BC_checks()) {
                this->zero_correction_related_to_equality_constrain(this->function(level), this->memory_.s[level]);
            }

            if (this->use_line_search() && (level == this->n_levels() - 1)) {
                auto ls_strategy_ = std::make_shared<utopia::Backtracking<Vector>>();
                Scalar alpha_ = 1.0;
                ls_strategy_->get_alpha(this->function(level),
                                        this->ml_derivs_.g[level],
                                        this->memory_.x[level],
                                        this->memory_.s[level],
                                        alpha_);
                std::cout << "alpha: " << alpha_ << "  \n";
                this->memory_.s[level] *= alpha_;
            }

            E_old = this->memory_.energy[level];
            this->memory_.x[level] += this->memory_.s[level];

            this->compute_s_global(level, this->memory_.s_working[level]);

            E_new = this->get_multilevel_energy(this->function(level), level, this->memory_.s_working[level]);

            //----------------------------------------------------------------------------
            //                        trial point acceptance
            //----------------------------------------------------------------------------
            // UTOPIA_NO_ALLOC_BEGIN("RMTR::region14");
            ared = E_old - E_new;
            rho = ared / coarse_reduction;
            if (!std::isfinite(E_new) || !std::isfinite(rho)) {
                rho = 0.0;
            }

            // in theory, this should never happen
            if (coarse_reduction <= 0) {
                rho = 0;
            }

            bool coarse_corr_taken = false;
            if (rho > this->rho_tol()) {
                coarse_corr_taken = true;
                this->memory_.energy[level] = E_new;

                // todo:: make sure that correct assumption
                this->get_multilevel_gradient(this->function(level), level, this->memory_.s_working[level], E_new);

                this->memory_.gnorm[level] = this->criticality_measure(level);
            } else {
                this->memory_.x[level] -= this->memory_.s[level];
                this->compute_s_global(level, this->memory_.s_working[level]);
            }

            //----------------------------------------------------------------------------
            //                                  trust region update
            //----------------------------------------------------------------------------
            converged = this->delta_update(rho, level, this->memory_.s_working[level]);
            if (converged && this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE) {
                std::cout << "cconverged after delta update  \n";
            }

            // because, x + Is_{l-1} does not have to be inside of the feasible set....
            // mostly case for rmtr_inf with bounds...
            if (rho > this->rho_tol() && converged == false) {
                converged = this->check_feasibility(level);
                if (this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && converged == true) {
                    std::cout << "- feasibility problems, terminating \n";
                }
            }

            if (this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && this->verbose() == true &&
                mpi_world_rank() == 0) {
                // just to see what is being printed
                std::string status = "RMTR_coarse_corr_stat, level: " + std::to_string(level);
                this->print_init_message(status,
                                         {" it. ",
                                          "   E_old     ",
                                          "   E_new",
                                          "ared   ",
                                          "  coarse_level_reduction  ",
                                          "  rho  ",
                                          "  delta ",
                                          "taken"});
                PrintInfo::print_iter_status(
                    this->_it_global,
                    {E_old, E_new, ared, coarse_reduction, rho, this->memory_.delta[level], Scalar(coarse_corr_taken)});
            }

            // UTOPIA_NO_ALLOC_END();

            // terminate, since TR rad. does not allow to take more corrections on given
            // level if (converged == true) {
            //     if (this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE) {
            //         std::cout << " converged last  \n";
            //     }
            //     return true;
            // }

            this->make_ml_iterate_feasible(level);

        } else if (mpi_world_rank() == 0 && this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE &&
                   this->verbose() == true) {
            std::cout << "--------- Recursion terminated due to non-smoothness of the "
                         "gradient, level: "
                      << level << " ----------------------- \n";
        }

        //----------------------------------------------------------------------------
        //                        postsmoothing
        //----------------------------------------------------------------------------

        if (this->post_smoothing_steps() != 0) {
            // auto post_smoothing_solve_type = (!smoothness_flg) ? COARSE_SOLVE :
            // POST_SMOOTHING;
            auto post_smoothing_solve_type = POST_SMOOTHING;
            this->local_tr_solve(level, post_smoothing_solve_type);
        }

        // exit(0);

        return true;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template <class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL>
    bool RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>::local_tr_solve(const SizeType &level,
                                                                     const LocalSolveType &solve_type) {
        SizeType it_success = 0, it = 0;
        Scalar ared = 0., pred = 0., rho = 0., energy_new = 9e9;
        bool make_grad_updates = true, make_hess_updates = true, converged = false, delta_converged = false;

        UTOPIA_NO_ALLOC_BEGIN("RMTR::region3");
        const bool exact_solve_flg = (solve_type == COARSE_SOLVE) ? true : false;

        this->initialize_local_solve(level, solve_type);

        make_hess_updates = this->init_deriv_loc_solve(this->function(level), level, solve_type);

        converged = this->check_local_convergence(it, it_success, level, this->memory_.delta[level], solve_type);

        if (this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && this->verbose() == true &&
            mpi_world_rank() == 0) {
            this->print_level_info(level);
            PrintInfo::print_iter_status(
                0,
                {this->memory_.gnorm[level], this->memory_.energy[level], ared, pred, rho, this->memory_.delta[level]});
        }

        UTOPIA_NO_ALLOC_END();

        it++;

        while (!converged) {
            UTOPIA_NO_ALLOC_BEGIN("RMTR::hessian_eval1");
            if (make_hess_updates) {
                bool update_h = (this->Hpost_lagging() && solve_type == POST_SMOOTHING) ? false : true;
                if (update_h) {
                    this->get_multilevel_hessian(this->function(level), level);
                }
            }
            UTOPIA_NO_ALLOC_END();

            //----------------------------------------------------------------------------
            //     solving constrained system to get correction and  building trial
            //     point
            //----------------------------------------------------------------------------
            // obtain correction
            UTOPIA_NO_ALLOC_BEGIN("RMTR::region5");
            this->solve_qp_subproblem(level, exact_solve_flg);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("RMTR::region6");
            // predicted reduction based on model
            pred = this->get_pred(level);

            // building trial point
            this->memory_.x[level] += this->memory_.s[level];

            this->compute_s_global(level, this->memory_.s_working[level]);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("RMTR::deriv_comp");
            energy_new = this->get_multilevel_energy(this->function(level), level, this->memory_.s_working[level]);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("RMTR::region701");
            ared = this->memory_.energy[level] - energy_new;

            rho = (ared < 0.0) ? 0.0 : ared / pred;
            rho = (rho != rho) ? 0.0 : rho;

            // just some safety checks
            if (!std::isfinite(energy_new) || !std::isfinite(rho) || has_nan_or_inf(this->memory_.x[level])) {
                rho = 0.0;
            }

            // update of hessian approx ...
            // TODO:: could be done in more elegant way....
            this->update_level(level, energy_new);

            //----------------------------------------------------------------------------
            //     acceptance of trial point
            //----------------------------------------------------------------------------
            // good reduction, accept trial point
            if (rho >= this->rho_tol()) {
                it_success++;
                this->memory_.energy[level] = energy_new;
                make_grad_updates = true;
            } else {
                this->memory_.x[level] -= this->memory_.s[level];  // return iterate into its initial state
                this->compute_s_global(level, this->memory_.s_working[level]);
                make_grad_updates = false;
            }

            //----------------------------------------------------------------------------
            //     updating level (deltas, hessian approx - new vectors, ...)
            //----------------------------------------------------------------------------
            delta_converged = this->delta_update(rho, level, this->memory_.s_working[level]);

            // TODO:: minimize norm computations
            // if(this->norm_schedule()==OUTER_CYCLE && this->verbosity_level() <
            // VERBOSITY_LEVEL_VERY_VERBOSE && (solve_type==POST_SMOOTHING || solve_type
            // == COARSE_SOLVE) && check_iter_convergence(it, it_success, level,
            // solve_type)) if(this->norm_schedule()==OUTER_CYCLE &&
            // this->verbosity_level() < VERBOSITY_LEVEL_VERY_VERBOSE &&
            // (solve_type==POST_SMOOTHING || solve_type == COARSE_SOLVE) &&
            // check_iter_convergence(it, it_success, level, solve_type))
            // {
            //     make_grad_updates = false;
            // }

            // can be more efficient, see commented lines below
            make_hess_updates = make_grad_updates;
            UTOPIA_NO_ALLOC_END();

            rho = this->update_local_grad(make_grad_updates, level, rho, energy_new);

            // else
            // {
            //     make_hess_updates = false;
            // }

            UTOPIA_NO_ALLOC_BEGIN("RMTR::region9");
            if (this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && this->verbose() == true &&
                mpi_world_rank() == 0) {
                PrintInfo::print_iter_status(it,
                                             {this->memory_.gnorm[level],
                                              this->memory_.energy[level],
                                              ared,
                                              pred,
                                              rho,
                                              this->memory_.delta[level]});
            }

            converged = (delta_converged == true) ? true
                                                  : this->check_local_convergence(
                                                        it, it_success, level, this->memory_.delta[level], solve_type);

            if (level == this->n_levels() - 1) {
                converged = (converged == true || this->memory_.gnorm[level] < this->atol()) ? true : false;
            }
            UTOPIA_NO_ALLOC_END();

            it++;
        }

        if (this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && this->verbose() == true &&
            mpi_world_rank() == 0) {
            std::cout << this->def_;
        }

        bool level_quit = ((this->criticality_measure_termination(level) == true) || delta_converged) ? true : false;

        return level_quit;
    }

}  // namespace utopia

#endif  // UTOPIA_RMTR_BASE_IMPL_HPP
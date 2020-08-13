#ifndef UTOPIA_DM_RMTR_SETUP_HPP
#define UTOPIA_DM_RMTR_SETUP_HPP

#include "utopia_BlockQPSolver.hpp"
#include "utopia_IPTransfer.hpp"
#include "utopia_Input.hpp"
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_MLIncrementalLoading.hpp"
#include "utopia_Multigrid.hpp"
#include "utopia_Multilevel.hpp"
#include "utopia_PFMassMatrix.hpp"
#include "utopia_RedundantQPSolver.hpp"
#include "utopia_make_unique.hpp"

#include <memory>

namespace utopia {

    // FIXME complete the overriding process
    template <class FunctionSpace, class ProblemType, class BCType, class ICType>
    class MLIncrementalLoading final : public IncrementalLoadingBase<FunctionSpace> {
    public:
        using Matrix = typename FunctionSpace::Matrix;
        using Vector = typename FunctionSpace::Vector;
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Super = utopia::IncrementalLoadingBase<FunctionSpace>;
        using Super::init;

        MLIncrementalLoading(FunctionSpace &space_coarse)
            : init_(false),
              n_levels_(2),
              n_coarse_sub_comm_(1),
              log_output_path_("rmtr_log_file.csv"),
              save_output_(true),
              mprgp_smoother_(false),
              hjsmn_smoother_(false) {
            spaces_.resize(2);
            spaces_[0] = make_ref(space_coarse);
        }

        void read(Input &in) override {
            IncrementalLoadingBase<FunctionSpace>::read(in);

            in.get("log_output_path", log_output_path_);
            in.get("n_coarse_sub_comm", n_coarse_sub_comm_);
            in.get("n_levels", n_levels_);
            in.get("save_output", save_output_);
            in.get("mprgp_smoother", mprgp_smoother_);
            in.get("hjsmn_smoother", hjsmn_smoother_);

            init_ml_setup();

            for (std::size_t l = 0; l < level_functions_.size(); l++) {
                level_functions_[l]->read(in);
                BC_conditions_[l]->read(in);
            }

            IC_->read(in);

            in.get("solver", *rmtr_);
            in.get("second_phase_ts", second_phase_time_stepper_);
        }

        bool init_ml_setup() {
            if (n_levels_ < 2) {
                std::cerr << "n_levels must be at least 2" << std::endl;
                return false;
            }

            spaces_.resize(n_levels_);

            level_functions_.resize(n_levels_);
            auto fun = std::make_shared<ProblemType>(*spaces_[0]);
            fun->use_crack_set_irreversibiblity(false);
            level_functions_[0] = fun;

            BC_conditions_.resize(n_levels_);
            auto bc = std::make_shared<BCType>(*spaces_[0]);
            BC_conditions_[0] = bc;

            transfers_.resize(n_levels_ - 1);

            for (SizeType i = 1; i < n_levels_; ++i) {
                spaces_[i] = spaces_[i - 1]->uniform_refine();

                auto fun = std::make_shared<ProblemType>(*spaces_[i]);

                if (i < n_levels_ - 1) {
                    fun->use_crack_set_irreversibiblity(false);
                } else {
                    fun->use_crack_set_irreversibiblity(true);
                }

                level_functions_[i] = fun;

                auto bc = std::make_shared<BCType>(*spaces_[i]);
                BC_conditions_[i] = bc;

                auto I = std::make_shared<Matrix>();
                spaces_[i - 1]->create_interpolation(*spaces_[i], *I);
                assert(!empty(*I));

                Matrix Iu;  // = *I;
                Iu.destroy();
                MatConvert(raw_type(*I), I->type_override(), MAT_INITIAL_MATRIX, &raw_type(Iu));
                Matrix R = transpose(Iu);

                PFMassMatrix<FunctionSpace> mass_matrix_assembler_fine(*spaces_[i]);
                Matrix M_fine;
                mass_matrix_assembler_fine.mass_matrix(M_fine);

                PFMassMatrix<FunctionSpace> mass_matrix_assembler_coarse(*spaces_[i - 1]);
                Matrix M_coarse;
                mass_matrix_assembler_coarse.mass_matrix(M_coarse);

                Matrix inv_lumped_mass = diag(1. / sum(M_coarse, 1));
                Matrix P = inv_lumped_mass * R * M_fine;

                transfers_[i - 1] = std::make_shared<IPTransferNested<Matrix, Vector>>(std::make_shared<Matrix>(Iu),
                                                                                       std::make_shared<Matrix>(P));
            }

            // initial conddition needs to be setup only on the finest level
            SizeType pf_comp = 0;
            IC_ = std::make_shared<ICType>(*spaces_.back(), pf_comp);

            //////////////////////////////////////////////// init solver
            ///////////////////////////////////////////////////
            // rmtr_ = std::make_shared<RMTR_inf<Matrix, Vector, TRKornhuberBoxKornhuber<Matrix, Vector>, SECOND_ORDER>
            // >(n_levels_); rmtr_ = std::make_shared<RMTR_inf<Matrix, Vector, TRBoundsGratton<Matrix, Vector>,
            // SECOND_ORDER> >(n_levels_);

            if (!rmtr_) {
                rmtr_ = std::make_shared<RMTR_inf<Matrix, Vector, TRBoundsGratton<Matrix, Vector>, SECOND_ORDER>>(
                    n_levels_);
            }

            // auto tr_strategy_fine   = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();
            // tr_strategy_fine->l1(true);

            std::shared_ptr<QPSolver<PetscMatrix, PetscVector>> tr_strategy_fine;

            if (mprgp_smoother_) {
                tr_strategy_fine = std::make_shared<utopia::MPGRP<Matrix, Vector>>();
            } else if (hjsmn_smoother_) {
                // auto qp = std::make_shared<SemismoothNewton<Matrix, Vector>>(
                //     std::make_shared<Factorization<Matrix, Vector>>());

                auto qp = std::make_shared<SemismoothNewton<Matrix, Vector>>(std::make_shared<MPGRP<Matrix, Vector>>());

                qp->max_it(2);
                // BlockQPSolver<Matrix, Vector> bqp(qp);
                tr_strategy_fine = std::make_shared<utopia::BlockQPSolver<Matrix, Vector>>(qp);
                // tr_strategy_fine->verbose(true);
            } else {
                tr_strategy_fine = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector>>();
            }

            std::shared_ptr<QPSolver<Matrix, Vector>> tr_strategy_coarse;
            if (n_coarse_sub_comm_ > 1 && n_coarse_sub_comm_ >= spaces_[0]->comm().size()) {
                spaces_[0]->comm().root_print("using redundant qp solver");
                auto qp = std::make_shared<utopia::MPGRP<Matrix, Vector>>();
                tr_strategy_coarse = std::make_shared<RedundantQPSolver<Matrix, Vector>>(qp, n_coarse_sub_comm_);
                // tr_strategy_coarse->verbose(true);
            } else {
                tr_strategy_coarse = std::make_shared<utopia::MPGRP<Matrix, Vector>>();
            }

            // auto ls = std::make_shared<GMRES<Matrix, Vector> >();
            // ls->pc_type("bjacobi");
            // auto tr_strategy_coarse = std::make_shared<utopia::SemismoothNewton<Matrix, Vector> >(ls);

            // rmtr_->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            rmtr_->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

            rmtr_->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr_->set_fine_tr_strategy(tr_strategy_fine);

            rmtr_->set_transfer_operators(transfers_);
            rmtr_->set_functions(level_functions_);
            rmtr_->verbose(true);

            init_ = true;

            return true;
        }

        FunctionSpace &fine_space() { return *spaces_.back(); }

        const FunctionSpace &fine_space() const { return *spaces_.back(); }

        std::shared_ptr<FunctionSpace> fine_space_ptr() { return spaces_.back(); }

        std::shared_ptr<const FunctionSpace> fine_space_ptr() const { return spaces_.back(); }

        ////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////

        void init_solution() override {
            spaces_.back()->create_vector(this->solution_);
            spaces_.back()->create_vector(this->lb_);
            rename("X", this->solution_);

            IC_->init(this->solution_);
            // if(this->use_pressure_ && !this->use_constant_pressure_){
            //     Vector & pressure_vec =  level_functions_.back()->pressure_field();
            //     // IC_->init(this->solution_, pressure_vec);
            // }

            spaces_.back()->apply_constraints(this->solution_);

            if (auto *fun_finest = dynamic_cast<ProblemType *>(level_functions_.back().get())) {
                fun_finest->old_solution(this->solution_);
            }

            // adding sol to all levels
            for (auto l = n_levels_ - 1; l > 0; l--) {
                auto *fun_fine = dynamic_cast<ProblemType *>(level_functions_[l].get());
                Vector &fine_sol = fun_fine->old_solution();

                auto *fun_coarse = dynamic_cast<ProblemType *>(level_functions_[l - 1].get());
                Vector &coarse_sol = fun_coarse->old_solution();
                spaces_[l]->create_vector(coarse_sol);

                transfers_[l - 1]->project_down(fine_sol, coarse_sol);

                std::cout<<"coarse_sol: "<< size(coarse_sol) << "  \n";

                // spaces_[l]->apply_constraints(coarse_sol);

                // transfers_[l]->restrict(fine_sol, coarse_sol);
            }
        }

        void write_to_file(FunctionSpace &space, const Scalar &time) override {
            if (save_output_) {
                // only finest level
                IncrementalLoadingBase<FunctionSpace>::write_to_file(space, time);

                // // all levels
                // for(auto l=0; l < spaces_.size(); l++){

                //     ProblemType * fun = dynamic_cast<ProblemType *>(level_functions_[l].get());
                //     Vector & sol  = fun->old_solution();
                //     rename("X", sol);

                //     std::cout<<"------ lev: "<< l << "   \n";
                //     std::cout<<"sol: "<< size(sol) << "  \n"; 

                //     auto name = this->output_path_+"_lev_"+ std::to_string(l)+"_time_"+std::to_string(time)+".vtk"; 
                //     std::cout<<"------ name: "<< name << " \n";

                //     spaces_[l]->write(name, sol);
                // }

                  // Vector sol = this->solution_; 


                  // for (auto l = n_levels_ - 1; l > 0; l--) 
                  // {

                  //     ProblemType * fun_coarse = dynamic_cast<ProblemType *>(level_functions_[l-1].get());
                  //     Vector sol_coarse  = 0*fun_coarse->get_eq_constrains_flg();

                  //     // transfers_[l - 1]->project_down(sol, sol_coarse);
                  //     transfers_[l - 1]->restrict(sol, sol_coarse);

                  //     // L2-fit projection 
                  //     // IPTransferNested<Matrix, Vector> * tr_nested = dynamic_cast<IPTransferNested<Matrix, Vector> *>(transfers_[l - 1].get());
                  //     // const Matrix & I = tr_nested->I(); 
                  //     // Vector s_help = transpose(I)*sol; 
                  //     // Matrix II = transpose(I) * I; 
                  //     // auto direct_solver = std::make_shared<LUDecomposition<PetscMatrix, PetscVector>>();
                  //     // direct_solver->solve(II, s_help, sol_coarse); 


                  //     spaces_[l-1]->apply_constraints(sol_coarse);
                  //     rename("X", sol_coarse);
                      

                  //     auto name = this->output_path_+"_lev_"+ std::to_string(l-1)+"_time_"+std::to_string(time)+".vtr"; 
                  //     spaces_[l-1]->write(name,
                  //     sol_coarse);


                  //     transfers_[l - 1]->interpolate(sol_coarse, sol); 
                  //     auto name2 = this->output_path_+"_lev_up_"+ std::to_string(l)+"_time_"+std::to_string(time)+".vtr"; 
                  //     rename("X", sol);
                  //     spaces_[l]->apply_constraints(sol);
                  //     spaces_[l]->write(name2, sol);


                  //     sol = sol_coarse; 
                  // }                

                Utopia::instance().set("log_output_path", log_output_path_);
            }
        }

        void prepare_for_solve() override {
            for (std::size_t l = 0; l < BC_conditions_.size(); l++) {
                BC_conditions_[l]->emplace_time_dependent_BC(this->time_);
            }

            // update fine level solution  and constraint
            spaces_.back()->apply_constraints(this->solution_);

            if (auto *fun_finest = dynamic_cast<ProblemType *>(level_functions_.back().get())) {
                fun_finest->build_irreversility_constraint(this->lb_);
            }

            for (std::size_t l = 0; l < BC_conditions_.size(); l++) {
                Vector &bc_flgs = level_functions_[l]->get_eq_constrains_flg();
                Vector &bc_values = level_functions_[l]->get_eq_constrains_values();

                spaces_[l]->apply_constraints(bc_values);
                spaces_[l]->build_constraints_markers(bc_flgs);

                // only on the finest level
                if (l == (BC_conditions_.size() - 1)) {
                    if (auto *fun_finest = dynamic_cast<ProblemType *>(level_functions_.back().get())) {
                        fun_finest->add_irr_values_markers(bc_values, bc_flgs);
                    }
                }

                // disp(bc_values, "bc_values");
                // disp(bc_flgs, "bc_flgs");
                level_functions_[l]->init_constraint_indices();
            }

            // if(this->use_pressure_){
            auto press_ts = this->pressure0_ + (this->time_ * this->pressure_increase_factor_);

            // if(this->use_constant_pressure_){
            // fe_problem_->setup_constant_pressure_field(press_ts);
            // std::cout<<"----- yes, constant pressure, "<< press_ts << " ......... \n";

            for (SizeType l = 0; l < n_levels_; l++) {
                auto *fun = dynamic_cast<ProblemType *>(level_functions_[l].get());
                fun->setup_constant_pressure_field(press_ts);
                fun->set_pressure(press_ts);
            }
            // }
            //     else{
            //         Vector & pressure_vec =  fe_problem_->pressure_field();
            //         // set_nonzero_elem_to(pressure_vec, press_ts);

            //         set_nonzero_elem_to(pressure_vec, (this->time_ * this->pressure_increase_factor_));
            // }
            // }
        }

        void update_time_step(const SizeType &conv_reason) override {
            if (this->adjust_dt_on_failure_ && conv_reason < 0) {
                if (auto *fun_finest = dynamic_cast<ProblemType *>(level_functions_.back().get())) {
                    fun_finest->get_old_solution(this->solution_);
                }

                // reset sol on all levels - important for BC conditions mostly s
                for (auto l = n_levels_ - 1; l > 0; l--) {
                    auto *fun_fine = dynamic_cast<ProblemType *>(level_functions_[l].get());
                    Vector &fine_sol = fun_fine->old_solution();

                    auto *fun_coarse = dynamic_cast<ProblemType *>(level_functions_[l - 1].get());
                    Vector &coarse_sol = fun_coarse->old_solution();

                    if (empty(coarse_sol)) {
                        spaces_[l]->create_vector(coarse_sol);
                    }

                    transfers_[l - 1]->project_down(fine_sol, coarse_sol);
                    spaces_[l]->apply_constraints(coarse_sol);
                    // transfers_[l]->restrict(fine_sol, coarse_sol);
                }

                this->time_ -= this->dt_;
                this->dt_ = this->dt_ * this->shrinking_factor_;
                this->time_ += this->dt_;
            } else {
                // std::cout<<"------- yes, updating...  \n";

                if (auto *fun_finest = dynamic_cast<ProblemType *>(level_functions_.back().get())) {
                    fun_finest->set_old_solution(this->solution_);
                }

                // update sol on all levels
                for (auto l = n_levels_ - 1; l > 0; l--) {
                    auto *fun_fine = dynamic_cast<ProblemType *>(level_functions_[l].get());
                    Vector &fine_sol = fun_fine->old_solution();

                    auto *fun_coarse = dynamic_cast<ProblemType *>(level_functions_[l - 1].get());
                    Vector &coarse_sol = fun_coarse->old_solution();

                    if (empty(coarse_sol)) {
                        spaces_[l]->create_vector(coarse_sol);
                    }

                    transfers_[l - 1]->project_down(fine_sol, coarse_sol);
                    spaces_[l]->apply_constraints(coarse_sol);
                    // transfers_[l]->restrict(fine_sol, coarse_sol);
                }

                if (this->pressure0_ != 0.0) {
                    this->write_to_file(*spaces_.back(), 1e-5 * this->time_);
                } else {
                    this->write_to_file(*spaces_.back(), this->time_);
                }

                if (this->time_ < second_phase_time_stepper_.start_time_) {
                    // increment time step
                    this->time_ += this->dt_;
                    this->time_step_counter_ += 1;
                } else {
                    this->time_ += second_phase_time_stepper_.dt_;
                    this->time_step_counter_ += 1;
                }
            }
        }

        void run() override {
            if (!init_) {
                init_ml_setup();
            }

            // init fine level spaces
            this->init(*spaces_[n_levels_ - 1]);

            this->time_step_counter_ = 0;
            while (this->time_ < this->final_time_) {
                if (mpi_world_rank() == 0) {
                    utopia::out() << "###################################################################### \n";
                    utopia::out() << "Time-step: " << this->time_step_counter_ << "  time:  " << this->time_
                                  << "  dt:  " << this->dt_ << " \n";
                    utopia::out() << "###################################################################### \n";
                }

                prepare_for_solve();

                // ////////////////////////////////////////////////////////////////////////////////////////////////////////
                rmtr_->solve(this->solution_);
                auto sol_status = rmtr_->solution_status();
                //////////////////////////////////////////////////////////////////////////////////////////////////////////
                update_time_step(sol_status.reason);
            }
        }

    private:
        bool init_;
        SizeType n_levels_;
        SizeType n_coarse_sub_comm_;

        std::vector<std::shared_ptr<FunctionSpace>> spaces_;
        std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> transfers_;

        std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector>>> level_functions_;
        std::vector<std::shared_ptr<BCType>> BC_conditions_;

        std::shared_ptr<ICType> IC_;
        std::string log_output_path_;

        std::shared_ptr<RMTR_inf<Matrix, Vector, TRBoundsGratton<Matrix, Vector>, SECOND_ORDER>> rmtr_;

        bool save_output_;

        TimeStepperInfo<Scalar> second_phase_time_stepper_;

        bool mprgp_smoother_;
        bool hjsmn_smoother_;
    };

}  // namespace utopia

#endif  // UTOPIA_DM_RMTR_SETUP_HPP

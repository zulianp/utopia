#ifndef UTOPIA_INCREMENTAL_LOADING_HPP
#define UTOPIA_INCREMENTAL_LOADING_HPP

#include "utopia_Base.hpp"
#include "utopia_RangeDevice.hpp"

// include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_BCSetup.hpp"
#include "utopia_Backtracking.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_InitialCondition.hpp"
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_LBFGS.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_MPRGP.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_PFMassMatrix.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_PseudoContinuationIncludes.hpp"
#include "utopia_QuasiNewtonBound.hpp"
#include "utopia_QuasiTrustRegionVariableBound.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_TrustRegionVariableBound.hpp"
// #include "utopia_petsc_Slepc.hpp"

#include "utopia_petsc_DM.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"

#include <chrono>
#include <cmath>
#include <random>

namespace utopia {

    template <class Scalar>
    class TimeStepperInfo : public Configurable {
    public:
        TimeStepperInfo() : dt_(1e-5), start_time_(9e9), final_time_(9e9) {}

        void read(Input &in) override {
            in.get("dt", dt_);
            in.get("start_time", start_time_);
            in.get("final_time", final_time_);
        }

    public:
        Scalar dt_;
        Scalar start_time_;
        Scalar final_time_;
    };

    template <class FunctionSpace>
    class IncrementalLoadingBase : public Configurable {
    public:
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        IncrementalLoadingBase()
            : dt_(1e-5),
              dt0_(1e-5),
              num_time_steps_(100),
              time_(0.0),
              final_time_(0.0),
              time_step_counter_(0.0),
              output_path_("out"),
              total_wall_clock_time_(0.0),

              shrinking_factor_(0.5),
              pressure0_(0.0),
              pressure_increase_factor_(0.0),
              increase_factor_(1.0),
              fracture_energy_time_stepping_(false)

        {}

        ~IncrementalLoadingBase() override = default;

        void read(Input &in) override {
            in.get("dt", dt_);
            in.get("num_time_steps", num_time_steps_);
            in.get("final_time", final_time_);
            in.get("output-path", output_path_);
            in.get("use_mprgp", use_mprgp_);
            in.get("shrinking_factor", shrinking_factor_);
            in.get("adjust_dt_on_failure", adjust_dt_on_failure_);
            in.get("pressure0", pressure0_);
            in.get("pressure_increase_factor", pressure_increase_factor_);
            in.get("use_pressure", use_pressure_);
            in.get("use_constant_pressure", use_constant_pressure_);
            in.get("fracture_energy_time_stepping", fracture_energy_time_stepping_);

            // E.P dynamic time stepping with fracture energy
            in.get("frac_energy_max_change", frac_energy_max_change_);
            in.get("frac_energy_min_change", frac_energy_min_change_);
            in.get("dt_min", dt_min_);
            in.get("dt_max", dt_max_);
            in.get("increase_factor", increase_factor_);

            // E.P Two timestepping
            in.get("use_two_time_steps", use_two_time_steps_);
            in.get("time_secondphase", time_secondphase_);
            in.get("dt_secondphase", dt_secondphase_);

            csv_file_name_ = this->output_path_ + "_energies.csv";
        }

        virtual void run() = 0;
        virtual void init_solution() = 0;
        virtual void prepare_for_solve() = 0;
        virtual void update_time_step(const SizeType &conv_reason, bool fracture_energy_monitoring) = 0;

        void init(FunctionSpace &space) {
            init_solution();

            if (pressure0_ != 0.0) {
                use_pressure_ = true;
            }

            write_to_file(space, 0.0);

            dt0_ = dt_;
            time_ = dt_;

            if (final_time_ == 0) {
                final_time_ = num_time_steps_ * dt_;
            }
        }

        // calculates global fracture energy, and compares it to the previous time step. Only returns true if max change
        // ,
        virtual bool must_reduce_time_step(Scalar &frac_energy) {
            if (this->frac_energy_old_ == 0.0) return false;  // dont check on the first time step

            if (frac_energy <= 0.0) return false;  // Do not reduce time step when no phase has propagated

            // Negative Old fracture energy treated as positive
            double frac_energy_old = this->frac_energy_old_ > 0.0 ? this->frac_energy_old_ : -this->frac_energy_old_;

            bool repeat = frac_energy / frac_energy_old > this->frac_energy_max_change_;

            if (mpi_world_rank() == 0) {
                utopia::out() << "-----------------------------------------------------"
                                 "----------------- \n";
                utopia::out() << "Fracture Energy -> New: " << frac_energy
                              << "\n                -> Old: " << this->frac_energy_old_ << "\nmeasured change ->  "
                              << frac_energy / frac_energy_old << "\nprescribed change-> "
                              << this->frac_energy_max_change_;
                utopia::out() << "\n-----------------------------------------------------"
                                 "----------------- \n";
            }

            // setting criteria for increasing time step ( a minimum increase in frac energy )
            this->increase_next_time_step_ =
                (frac_energy / frac_energy_old < this->frac_energy_min_change_ && frac_energy / frac_energy_old > 0.95);

            return repeat;
        }

    protected:
        virtual void write_to_file(FunctionSpace &space, const Scalar &time) {
            space.write(output_path_ + "_" + std::to_string(time) + ".vtr", solution_);
        }

    protected:
        Scalar dt_, dt0_;
        Scalar num_time_steps_;
        Scalar time_;
        Scalar final_time_;
        SizeType time_step_counter_;
        std::string output_path_;
        bool use_mprgp_{false};
        Scalar total_wall_clock_time_;

        bool adjust_dt_on_failure_{true};
        Scalar shrinking_factor_;
        Scalar pressure0_;
        Scalar pressure_increase_factor_;
        bool use_pressure_{true};
        bool use_constant_pressure_{false};

        Scalar frac_energy_old_{0.0};
        Scalar frac_energy_max_change_{1e10};
        Scalar dt_min_{0}, dt_max_{1e3};
        Scalar frac_energy_min_change_{1e10};
        bool increase_next_time_step_{false};

        bool use_two_time_steps_{false};
        Scalar dt_secondphase_{1};
        Scalar time_secondphase_{1};

        Scalar increase_factor_{1};

        Vector solution_;
        Vector lb_;  // this is quite particular for PF-frac        E.P Question: What is this?
        Vector ub_;  // this is quite particular for PF-frac

        // E.P For Exporting energy file and iterations
        std::string csv_file_name_;
        bool fracture_energy_time_stepping_ = false;
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template <class FunctionSpace, class ProblemType>
    class IncrementalLoading : public IncrementalLoadingBase<FunctionSpace> {
    public:
        // expose inner types
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;
        using Matrix = typename FunctionSpace::Matrix;

        IncrementalLoading(FunctionSpace &space, InitialCondition<FunctionSpace> &IC, BCSetup<FunctionSpace> &BC)
            : space_(space), IC_(IC), BC_(BC) {
            fe_problem_ = std::make_shared<ProblemType>(space_);
        }

        void read(Input &in) override {
            IncrementalLoadingBase<FunctionSpace>::read(in);
            IC_.read(in);
            BC_.read(in);

            fe_problem_->read(in);

            in.get("use_box_constraints", use_box_constraints_);

            choose_solver();

            if (use_box_constraints_ == false) in.get("solver", *tr_solver_);

            if (use_box_constraints_ == true) in.get("solver", *tr_solver_box_);
        }

        void choose_solver() {
            // if (tr_solver_) return;
            // if (tr_solver_box_) return;

            UTOPIA_TRACE_REGION_BEGIN("IncrementalLoading::choose_solver(...)");

            std::cout << "incremental loading ..... \n";

            // exit(0);

            if (use_box_constraints_) {
                std::shared_ptr<QPSolver<PetscMatrix, PetscVector>> qp_solver;

                if (this->use_mprgp_) {
                    // MPRGP sucks as a solver, as it can not be preconditioned easily ...
                    qp_solver = std::make_shared<utopia::MPRGP<PetscMatrix, PetscVector>>();
                } else {
                    // tao seems to be faster until it stalls ...
                    // auto linear_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
                    // auto linear_solver = std::make_shared<SteihaugToint<PetscMatrix, PetscVector>>();
                    auto linear_solver = std::make_shared<Lanczos<PetscMatrix, PetscVector>>();
                    linear_solver->pc_type("bjacobi");
                    qp_solver = std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector>>(linear_solver);
                }

                tr_solver_box_ = std::make_shared<TrustRegionVariableBound<PetscMatrix, PetscVector>>(qp_solver);
            } else {
#if 0
                    auto qp_solver = std::make_shared<utopia::SteihaugToint<PetscMatrix, PetscVector, HOMEMADE>>();
                    auto prec = std::make_shared<KSPSolver<PetscMatrix, PetscVector>>();
                    prec->ksp_type("bcstab");
                    prec->pc_type("hypre");
                    qp_solver->set_preconditioner(prec);
#else
                auto qp_solver = std::make_shared<utopia::Lanczos<PetscMatrix, PetscVector>>();
                qp_solver->pc_type("hypre");
#endif

                tr_solver_ = std::make_shared<TrustRegion<PetscMatrix, PetscVector>>(qp_solver);
            }

            UTOPIA_TRACE_REGION_END("IncrementalLoading::choose_solver(...)");
        }

        void init_solution() override {
            UTOPIA_TRACE_REGION_BEGIN("IncrementalLoading::init_solution(...)");

            space_.create_vector(this->solution_);
            space_.create_vector(this->lb_);
            rename("X", this->solution_);

            IC_.init(this->solution_);
            if (this->use_pressure_ && !this->use_constant_pressure_) {
                Vector &pressure_vec = fe_problem_->pressure_field();
                IC_.init(this->solution_, pressure_vec);
            }

            space_.apply_constraints(this->solution_);
            fe_problem_->set_old_solution(this->solution_);
            fe_problem_->write_to_file(
                this->output_path_, this->solution_, this->time_);  // outputting initial stress file

            UTOPIA_TRACE_REGION_END("IncrementalLoading::init_solution(...)");
        }

        void prepare_for_solve() override {
            UTOPIA_TRACE_REGION_BEGIN("IncrementalLoading::prepare_for_solve(...)");

            BC_.emplace_time_dependent_BC(this->time_);
            space_.apply_constraints(this->solution_);
            //fe_problem_->set_dt(this->dt_); E.P Removed, was only needed for mobility
            fe_problem_->set_time(this->time_);

            if (this->use_pressure_) {
                auto press_ts = this->pressure0_ + (this->time_ * this->pressure_increase_factor_);

                if (this->use_constant_pressure_) {
                    fe_problem_->setup_constant_pressure_field(press_ts);
                } else {
                    Vector &pressure_vec = fe_problem_->pressure_field();
                    // set_nonzero_elem_to(pressure_vec, press_ts);

                    set_nonzero_elem_to(pressure_vec, (this->time_ * this->pressure_increase_factor_));
                }
            }

            // E.P CHECK - What is this doing here.
            fe_problem_->build_irreversility_constraint(this->lb_);

            UTOPIA_TRACE_REGION_END("IncrementalLoading::prepare_for_solve(...)");
        }

        void update_time_step(const SizeType &conv_reason, bool fracture_energy_monitoring = false) override {
            UTOPIA_TRACE_REGION_BEGIN("IncrementalLoading::update_time_step(...)");

            if (!fracture_energy_monitoring) {
                if (this->adjust_dt_on_failure_ && conv_reason < 0) {
                    // reset solution
                    fe_problem_->get_old_solution(this->solution_);

                    this->time_ -= this->dt_;
                    this->dt_ = this->dt_ * this->shrinking_factor_;
                    this->time_ += this->dt_;
                } else if (this->use_two_time_steps_) {

                    //update solution
                    fe_problem_->set_old_solution(this->solution_);

                    // Writing strain stress solution to file
                    fe_problem_->write_to_file(this->output_path_, this->solution_, this->time_);

                    // Modify future time step
                    if (this->time_ >= this->time_secondphase_)
                        this->dt_ = this->dt_secondphase_;

                    //prepare new solution
                    this->time_ += this->dt_;
                    this->time_step_counter_ += 1;


                } else {
                    rename("X", this->solution_);

                    if (this->pressure0_ != 0.0) {
                        this->write_to_file(space_, 1e-5 * this->time_);
                    } else {
                        this->write_to_file(space_, this->time_);
                    }

                    fe_problem_->set_old_solution(this->solution_);

                    // Writing strain stress solution to file
                    fe_problem_->write_to_file(this->output_path_, this->solution_, this->time_);

                    // increment time step
                    this->time_ += this->dt_;
                    this->time_step_counter_ += 1;
                }
            } else {
                Scalar trial_fracture_energy{0.0};
                bool repeat_step = this->adjust_dt_on_failure_ && conv_reason < 0;
                if (!repeat_step) {
                    // E.P Calc frac energy at new solution
                    // And repeat if larger than a certain tolerance
                    fe_problem_->fracture_energy(this->solution_, trial_fracture_energy);
                    repeat_step = this->must_reduce_time_step(
                        trial_fracture_energy);  // this also checks if frac energy change is too small

                    // if we are at the minimum time step, do not shrink
                    if (repeat_step == true && this->dt_ * this->shrinking_factor_ < this->dt_min_) {
                        repeat_step = false;
                        if (mpi_world_rank() == 0) {
                            utopia::out() << "time step " << this->dt_ << " too close to minumum " << this->dt_min_
                                          << " not reducing\n";
                        }
                    }
                }  // repeat step calibrated

                if (repeat_step) {  // repeat time step
                    if (mpi_world_rank() == 0) {
                        utopia::out() << "------- Repeating time step\n";
                    }
                    this->time_ -= this->dt_;
                    this->dt_ = this->dt_ * this->shrinking_factor_;
                    this->time_ += this->dt_;

                    fe_problem_->get_old_solution(this->solution_);
                    //fe_problem_->set_dt(this->dt_); E.P Removed, was only needed for mobility

                } else {  // Advance time step

                    fe_problem_->set_old_solution(this->solution_);
                    //fe_problem_->set_dt(this->dt_); E.P Removed, was only needed for mobility

                    if (this->pressure0_ != 0.0) {
                        this->write_to_file(space_, 1e-5 * this->time_);
                    }

                    // Advance time step
                    this->time_step_counter_ += 1;

                    // Choosing time step
                    if (this->increase_next_time_step_ && this->dt_ * this->increase_factor_ < this->dt_max_) {
                        this->dt_ *= this->increase_factor_;     // E.P Increasing time step if decided earlier by
                                                                 // IncrementalLoadingBase
                        this->increase_next_time_step_ = false;  // reset to zero for next time step
                        if (mpi_world_rank() == 0) {
                            utopia::out() << "------- Increasing next time step:  " << this->dt_ << "\n";
                        }
                    }

                    // E.P this also updates the fracture energy at old time step
                    // Need to call for dynamic time stepping strategy that uses fracture energy
                    this->frac_energy_old_ = trial_fracture_energy;  // only store old value if we dont repeat

                    // Writing strain stress solution to file
                    fe_problem_->write_to_file(this->output_path_, this->solution_, this->time_);
                    export_energies_csv(trial_fracture_energy);

                    // Advancing time step (no Second Phase - E.P Taken away )
                    this->time_ += this->dt_;  // increment time step
                }
            }

            UTOPIA_TRACE_REGION_END("IncrementalLoading::update_time_step(...)");
        }

        void export_energies_csv(Scalar fracture_energy) {
            if (!this->csv_file_name_.empty()) {
                CSVWriter writer{};
                Scalar elastic_energy = 0.0, ela_en_mid = 0.0, fra_en_mid = 0.0, tcv = 0.0, error_cod = 0.0,
                       residual = 0.0, iterations = 0.0;

                fe_problem_->elastic_energy(this->solution_, elastic_energy);
                fe_problem_->elastic_energy_in_middle_layer(this->solution_, ela_en_mid);
                //                fe_problem_->fracture_energy_in_middle_layer(this->solution_, fra_en_mid);

                fe_problem_->compute_tcv(this->solution_, tcv);

                if (!use_box_constraints_) {
                    residual = tr_solver_->get_gnorm();
                    iterations = tr_solver_->get_iterations();
                }

                if (FunctionSpace::Dim == 3) {
                    fe_problem_->compute_cod(this->solution_, error_cod);
                }

                if (mpi_world_rank() == 0) {
                    if (!writer.file_exists(this->csv_file_name_)) {
                        writer.open_file(this->csv_file_name_);
                        writer.write_table_row<std::string>({"time step",
                                                             "time",
                                                             "elastic_energy",
                                                             "fracture_energy",
                                                             "elast_en_mid_layer",
                                                             "frac_en_mid_layer",
                                                             "total_crack_vol",
                                                             "gnorm",
                                                             "iterations"});
                    } else {
                        writer.open_file(this->csv_file_name_);
                    }

                    writer.write_table_row<Scalar>({Scalar(this->time_step_counter_),
                                                    this->time_,
                                                    elastic_energy,
                                                    fracture_energy,
                                                    ela_en_mid,
                                                    fra_en_mid,
                                                    tcv,
                                                    residual,
                                                    iterations});
                    writer.close_file();
                }
            }
        }

        // allow passing solver
        void run() override {
            UTOPIA_TRACE_REGION_BEGIN("IncrementalLoading::run(...)");

            // just for testing purposes...
            // fe_problem_->use_crack_set_irreversibiblity(true);
            //            fe_problem_->turn_off_cu_coupling(true);
            //            fe_problem_->turn_off_uc_coupling(true);

            this->init(space_);

            while (this->time_ < this->final_time_) {
                if (mpi_world_rank() == 0) {
                    utopia::out() << "###################################################################### \n";
                    utopia::out() << "Time-step: " << this->time_step_counter_ << "  time:  " << this->time_
                                  << "  dt:  " << this->dt_ << " \n";
                    utopia::out() << "###################################################################### \n";
                    utopia::out() << "use_box_constraints_  " << use_box_constraints_ << "  \n";
                }

                // fe problem is missing
                prepare_for_solve();

                if (this->time_step_counter_ == 1) {
                    fe_problem_->set_old_solution(this->solution_);
                    //fe_problem_->set_dt(this->dt_); //E.P removed - was only for used for mobility
                    fe_problem_->export_material_params(this->output_path_);
                }

                ////////////////////////////////////////////////////////////////////////////////////////////////////////
                SizeType conv_reason;

                if (use_box_constraints_) {
                    // auto box = make_lower_bound_constraints(make_ref(this->lb_));

                    PetscVector ub = 0.0 * this->lb_;
                    ub.set(1.0);

                    auto box = make_box_constaints(make_ref(this->lb_), make_ref(ub));
                    tr_solver_box_->set_box_constraints(box);
                    // tr_solver_->atol(1e-6);
                    // tr_solver_->max_it(200);

                    tr_solver_box_->solve(*fe_problem_, this->solution_);

                    auto sol_status = tr_solver_box_->solution_status();
                    conv_reason = sol_status.reason;
                } else {
                    // disp(this->solution_, "this->solution_");

                    tr_solver_->solve(*fe_problem_, this->solution_);
                    auto sol_status = tr_solver_->solution_status();
                    conv_reason = sol_status.reason;
                }

                // ////////////////////////////////////////////////////////////////////////////////////////////////////////

                // auto hess_approx = std::make_shared<JFNK<utopia::PetscVector>>(*fe_problem_);
                // auto qp_solver = std::make_shared<utopia::MPGRP<PetscMatrix, PetscVector>>();
                // qp_solver->max_it(1000);

                // QuasiTrustRegionVariableBound<utopia::PetscVector> tr_solver(hess_approx, qp_solver);

                // // tr_solver.max_it(10);
                // tr_solver.verbose(true);

                // PetscVector ub = 0.0 * this->lb_;
                // ub.set(1.0);
                // auto box = make_box_constaints(make_ref(this->lb_), make_ref(ub));
                // tr_solver.set_box_constraints(box);

                // tr_solver.solve(*fe_problem_, this->solution_);
                // auto sol_status = tr_solver.solution_status();
                // const auto conv_reason = sol_status.reason;

                //////////////////////////////////////////////////////////////////////////////////////////////////////////
                // auto linear_solver = std::make_shared<GMRES<PetscMatrix, PetscVector>>();
                // linear_solver->atol(1e-14);
                // linear_solver->max_it(10000);
                // linear_solver->pc_type(PCILU);

                // // AffineSimilarity<PetscMatrix, PetscVector> solver(linear_solver);
                // ASTRUM<PetscMatrix, PetscVector> solver(linear_solver);
                // // PseudoContinuation<PetscMatrix, PetscVector> solver(linear_solver);

                // // assemble mass matrix
                // PFMassMatrix<FunctionSpace> mass_matrix_assembler(space_);
                // PetscMatrix M;
                // mass_matrix_assembler.mass_matrix(M);
                // solver.set_mass_matrix(M);

                // solver.verbose(true);
                // solver.atol(1e-5);
                // solver.stol(1e-12);
                // // solver.max_it(300);

                // PetscVector ub = 0.0 * this->lb_;
                // ub.set(1.0);
                // auto box = make_box_constaints(make_ref(this->lb_), make_ref(ub));
                // solver.set_box_constraints(box);

                // solver.solve(*fe_problem_, this->solution_);
                // const auto conv_reason = 1;
                ////////////////////////////////////////////////////////////////////////////////////////////////////////

                update_time_step(conv_reason, this->fracture_energy_time_stepping_);
            }

            UTOPIA_TRACE_REGION_END("IncrementalLoading::run(...)");
        }

    public:  // made public because of nvcc
        void set_nonzero_elem_to(Vector &v, const Scalar &val) { v.shift(val); }

    private:
        FunctionSpace &space_;
        InitialCondition<FunctionSpace> &IC_;
        BCSetup<FunctionSpace> &BC_;

        bool use_box_constraints_{false};

        std::shared_ptr<ProblemType> fe_problem_;
        std::shared_ptr<TrustRegionVariableBound<Matrix, Vector>> tr_solver_box_;
        std::shared_ptr<TrustRegion<Matrix, Vector>> tr_solver_;
    };

}  // namespace utopia

#endif

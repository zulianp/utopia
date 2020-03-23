#ifndef UTOPIA_INCREMENTAL_LOADING_HPP
#define UTOPIA_INCREMENTAL_LOADING_HPP


#include "utopia_Base.hpp"
#include "utopia_RangeDevice.hpp"

//include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_Core.hpp"
#include "utopia_PetscDM.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_petsc.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_petsc_dma_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_MPRGP.hpp"
#include "utopia_TrustRegionVariableBound.hpp"
#include "utopia_QuasiNewtonBound.hpp"
#include "utopia_Backtracking.hpp"
#include "utopia_LBFGS.hpp"
#include "utopia_QuasiTrustRegionVariableBound.hpp"
#include "utopia_InitialCondition.hpp"
#include "utopia_BCSetup.hpp"


#include <random>
#include <cmath>
#include <chrono>

namespace utopia {

    template<class FunctionSpace>
    class IncrementalLoadingBase : public Configurable
    {
        public:
            using SizeType       = typename FunctionSpace::SizeType;
            using Scalar         = typename FunctionSpace::Scalar;
            using Vector         = typename FunctionSpace::Vector;

            IncrementalLoadingBase():
            dt_(1e-5),
            dt0_(1e-5),
            num_time_steps_(100),
            time_(0.0),
            output_path_("out.vtr"), 
            adjust_dt_on_failure_(true), 
            shrinking_factor_(0.5), 
            pressure0_(0.0),
            pressure_increase_factor_(0.0),
            use_pressure_(false),
            use_constant_pressure_(false)
            {

            }


            virtual ~IncrementalLoadingBase()
            {

            }


            virtual void read(Input &in) override
            {
                in.get("dt", dt_);
                in.get("num_time_steps", num_time_steps_);
                in.get("output-path", output_path_);
                in.get("use_mprgp", use_mprgp_);
                in.get("shrinking_factor", shrinking_factor_);
                in.get("adjust_dt_on_failure", adjust_dt_on_failure_);
                in.get("pressure0", pressure0_);
                in.get("pressure_increase_factor", pressure_increase_factor_);
                in.get("use_pressure", use_pressure_);
                in.get("use_constant_pressure", use_constant_pressure_);
            }

            // allow passing solver
            virtual void run(Input &in) = 0;
            virtual void init_solution() = 0;
            virtual void prepare_for_solve() = 0;
            virtual void update_time_step(const SizeType & conv_reason) = 0;


            void init(Input &in, FunctionSpace & space)
            {
                read(in);
                init_solution();

                if(pressure0_ != 0.0){
                    use_pressure_ = true;
                }


                write_to_file(space, 0.0);

                dt0_ = dt_;
                time_ = dt_;
            }


        protected:
            virtual void write_to_file(FunctionSpace & space, const Scalar & time)
            {
                space.write(output_path_+"_"+std::to_string(time)+".vtr", solution_);
            }



        protected:
            Scalar dt_, dt0_;
            Scalar num_time_steps_;
            Scalar time_;
            std::string output_path_;
            bool use_mprgp_;
            bool adjust_dt_on_failure_;
            Scalar shrinking_factor_;
            Scalar pressure0_;
            Scalar pressure_increase_factor_;
            bool use_pressure_;
            bool use_constant_pressure_;

            Vector solution_;
            Vector lb_; // this is quite particular for PF-frac
    };





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template<class FunctionSpace, class ProblemType>
    class IncrementalLoading : public IncrementalLoadingBase<FunctionSpace>
    {
        public:
            //expose inner types
            using SizeType       = typename FunctionSpace::SizeType;
            using Scalar         = typename FunctionSpace::Scalar;
            using Vector         = typename FunctionSpace::Vector;

            IncrementalLoading(FunctionSpace & space, InitialCondition<FunctionSpace> & IC, BCSetup<FunctionSpace> & BC):
            space_(space),
            IC_(IC),
            BC_(BC)
            {

            }


            void read(Input &in) override
            {
                IncrementalLoadingBase<FunctionSpace>::read(in);
                IC_.read(in);
                BC_.read(in);
                fe_problem_->read(in);
            }

            void init_solution() override
            {
                space_.create_vector(this->solution_);
                space_.create_vector(this->lb_);
                rename("X", this->solution_);

                IC_.init(this->solution_);
                if(this->use_pressure_ && !this->use_constant_pressure_){
                    Vector & pressure_vec =  fe_problem_->pressure_field();
                    IC_.init(this->solution_, pressure_vec);
                }

                space_.apply_constraints(this->solution_);
                fe_problem_->old_solution(this->solution_);
            }


            void prepare_for_solve() override
            {
                BC_.emplace_time_dependent_BC(this->time_);
                space_.apply_constraints(this->solution_);

                if(this->use_pressure_){
                    auto press_ts = this->pressure0_ + (this->time_ * this->pressure_increase_factor_);

                    if(this->use_constant_pressure_){
                        fe_problem_->setup_constant_pressure_field(press_ts);
                    }
                    else{
                        Vector & pressure_vec =  fe_problem_->pressure_field();
                        // set_nonzero_elem_to(pressure_vec, press_ts);

                        set_nonzero_elem_to(pressure_vec, (this->time_ * this->pressure_increase_factor_));
                    }
                }

                fe_problem_->build_irreversility_constraint(this->lb_);
            }



            void update_time_step(const SizeType & conv_reason) override
            {
                if(this->adjust_dt_on_failure_ && conv_reason < 0){
                        // reset solution
                        fe_problem_->get_old_solution(this->solution_);

                        this->time_ -= this->dt_;
                        this->dt_ = this->dt_ * this->shrinking_factor_;
                        this->time_ += this->dt_;
                    }
                    else{
                        rename("X", this->solution_);

                        if(this->pressure0_!= 0.0){
                            this->write_to_file(space_, 1e-5*this->time_);
                        }
                        else{
                            this->write_to_file(space_, this->time_);
                        }

                        fe_problem_->old_solution(this->solution_);

                        // increment time step
                        this->time_ += this->dt_;
                    }
            }



            // allow passing solver
            void run(Input &in) override
            {
                fe_problem_ = std::make_shared<ProblemType>(space_);
                this->init(in, space_);


                for (auto t=1; t < this->num_time_steps_; t++)
                {
                    if(mpi_world_rank()==0){
                        std::cout<<"###################################################################### \n";
                        std::cout<<"Time-step: "<< t << "  time:  "<< this->time_ << "  dt:  "<< this->dt_ << " \n";
                        std::cout<<"###################################################################### \n";
                    }

                    // fe problem is missing
                    prepare_for_solve();


                    ////////////////////////////////////////////////////////////////////////////////////////////////////////
                    std::shared_ptr<QPSolver<PetscMatrix, PetscVector>> qp_solver;
                    if(this->use_mprgp_) {
                        // MPRGP sucks as a solver, as it can not be preconditioned easily ...
                        qp_solver = std::make_shared<utopia::MPGRP<PetscMatrix, PetscVector> >();
                    } else {
                        // tao seems to be faster until it stalls ...
                        // auto linear_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
                        auto linear_solver = std::make_shared<GMRES<PetscMatrix, PetscVector>>();
                        linear_solver->max_it(200);
                        linear_solver->pc_type("bjacobi");
                        qp_solver = std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector> >(linear_solver);
                    }

                    qp_solver->max_it(1000);
                    TrustRegionVariableBound<PetscMatrix, PetscVector> solver(qp_solver);
                    auto box = make_lower_bound_constraints(make_ref(this->lb_));
                    solver.set_box_constraints(box);
                    in.get("solver", solver);
                    solver.solve(*fe_problem_, this->solution_);
                    auto sol_status = solver.solution_status();

                    ////////////////////////////////////////////////////////////////////////////////////////////////////////

                    const auto conv_reason = sol_status.reason;
                    update_time_step(conv_reason);

                }

            }



    public: //made public because of nvcc

        void set_nonzero_elem_to(Vector & v, const Scalar & val)
        {
            parallel_transform(v, UTOPIA_LAMBDA(const SizeType &i, const Scalar &vi) -> Scalar
            {
                return vi + val;
            });
        }



        private:
            FunctionSpace & space_;
            InitialCondition<FunctionSpace> & IC_;
            BCSetup<FunctionSpace> & BC_;

            std::shared_ptr<ProblemType > fe_problem_;
    };



}


#endif
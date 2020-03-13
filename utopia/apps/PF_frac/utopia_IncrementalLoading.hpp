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


#include <random>
#include <cmath>
#include <chrono>

namespace utopia {

    template<class FunctionSpace>
    class BCSetup : public Configurable
    {
        public:
            using Scalar    = typename FunctionSpace::Scalar;

            BCSetup(FunctionSpace & space):
            space_(space)
            {

            }

            virtual ~BCSetup()
            {

            }

            void read(Input &in) override
            {

            }            

            virtual void emplace_time_dependent_BC(const Scalar & time) = 0; 

        protected:
            FunctionSpace & space_; 
    };


    template<class FunctionSpace>
    class PFFracFixAllDisp2D : public BCSetup<FunctionSpace>
    {
        public:
            using Scalar    = typename FunctionSpace::Scalar;
            using Vector   = typename FunctionSpace::Vector;

            PFFracFixAllDisp2D(FunctionSpace & space): BCSetup<FunctionSpace>(space)
            {

            }

            void emplace_time_dependent_BC(const Scalar & /*time*/) override
            {
                static const int Dim   = FunctionSpace::Dim;

                using Point          = typename FunctionSpace::Point;
                this->space_.reset_bc(); 

                for(int d = 1; d < Dim + 1; ++d) {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::left(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                        return 0.0;
                        },
                        d
                    );

                    this->space_.emplace_dirichlet_condition(
                        SideSet::right(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                            return 0.0;
                        },
                        d
                        );

                    this->space_.emplace_dirichlet_condition(
                        SideSet::top(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                            return 0.0;
                        },
                        d
                        ); 

                    this->space_.emplace_dirichlet_condition(
                        SideSet::bottom(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                            return 0.0;
                        },
                        d
                        );
                }
            }
    };   



    template<class FunctionSpace>
    class PFFracTension2D : public BCSetup<FunctionSpace>
    {
        public:
            using Scalar    = typename FunctionSpace::Scalar;
            using Vector    = typename FunctionSpace::Vector;

            PFFracTension2D(FunctionSpace & space, const Scalar & disp_y=1.0): BCSetup<FunctionSpace>(space), disp_y_(disp_y)
            {

            }

            void read(Input &in) override
            {
                in.get("disp_y", disp_y_);   
            }               

            void emplace_time_dependent_BC(const Scalar & time) override
            {
                static const int Dim   = FunctionSpace::Dim;

                using Point          = typename FunctionSpace::Point;
                this->space_.reset_bc(); 

                this->space_.emplace_dirichlet_condition(
                    SideSet::left(),
                    UTOPIA_LAMBDA(const Point &p) -> Scalar {
                        return 0.0;
                    },
                    1
                    );

                this->space_.emplace_dirichlet_condition(
                    SideSet::right(),
                    UTOPIA_LAMBDA(const Point &p) -> Scalar {
                        return disp_y_ * time;
                    },
                    1
                    );

                for(int d = 2; d < Dim + 1; ++d) {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::left(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                        return 0.0;
                        },
                        d
                    );

                    this->space_.emplace_dirichlet_condition(
                        SideSet::right(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                            return 0.0;
                        },
                        d
                        );
                    }
            }

            private:
                Scalar disp_y_; 
    };        




    template<class FunctionSpace>
    class IncrementalLoading : public Configurable
    {
        public: 

            static const int Dim   = FunctionSpace::Dim;
            static const int NVars = FunctionSpace::Dim + 1;

            //expose inner types
            using Mesh           = typename FunctionSpace::Mesh;
            using Elem           = typename FunctionSpace::Shape;
            using ElemView       = typename FunctionSpace::ViewDevice::Elem;
            using SizeType       = typename FunctionSpace::SizeType;
            using Scalar         = typename FunctionSpace::Scalar;
            using Dev            = typename FunctionSpace::Device;
            using Point          = typename FunctionSpace::Point;
            using Vector         = typename FunctionSpace::Vector;
            using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;

            static const int NNodes = Elem::NNodes;

            using FEFunction     = utopia::FEFunction<FunctionSpace>;
            using Quadrature     = utopia::Quadrature<Elem, 2>;
            using Parameters     = typename IsotropicPhaseFieldForBrittleFractures<FunctionSpace>::Parameters;

            IncrementalLoading(FunctionSpace & space, InitialCondition<FunctionSpace> & IC, BCSetup<FunctionSpace> & BC): 
            space_(space), 
            IC_(IC), 
            BC_(BC), 
            dt_(1e-5), 
            dt0_(1e-5),
            num_time_steps_(100), 
            time_(0.0),
            output_path_("out.vtk"), 
            use_mprgp_(false), 
            adjust_dt_on_failure_(true), 
            shrinking_factor_(0.5), 
            pressure0_(0.0),
            pressure_increase_factor_(0.0)
            {

            }


            void read(Input &in) override
            {
                in.get("dt", dt_);   
                in.get("num_time_steps", num_time_steps_);   
                in.get("output-path", output_path_);
                in.get("use_mprgp", use_mprgp_); 
                in.get("shrinking_factor", shrinking_factor_); 
                in.get("adjust_dt_on_failure", adjust_dt_on_failure_); 
                in.get("pressure0", pressure0_); 
                in.get("pressure_increase_factor", pressure_increase_factor_);
            }

            // allow passing solver 
            template <class ProblemType>
            void run(Input &in)
            {
                this->read(in);  
                IC_.read(in); 
                BC_.read(in); 

                ProblemType fe_problem(space_); 
                fe_problem.read(in);  

                space_.create_vector(solution_);
                space_.create_vector(lb_);
                rename("X", solution_);

                IC_.init(solution_); 

                space_.apply_constraints(solution_);                    
                fe_problem.old_solution(solution_); 

                space_.write(output_path_+"_"+std::to_string(0.0)+".vtk", solution_);     

                dt0_ = dt_; 
                time_ = dt_; 

                for (auto t=1; t < num_time_steps_; t++)
                {
                    if(mpi_world_rank()==0){
                        std::cout<<"###################################################################### \n"; 
                        std::cout<<"Time-step: "<< t << "  time:  "<< time_ << "  dt:  "<< dt_ << " \n"; 
                        std::cout<<"###################################################################### \n"; 
                    }
             
                    BC_.emplace_time_dependent_BC(time_); 
                    space_.apply_constraints(solution_);    

                    if(pressure0_!= 0.0){
                        fe_problem.set_pressure(pressure0_ + (time_ * pressure_increase_factor_) );     
                    }
                    
                    fe_problem.build_irreversility_constraint(lb_); 

                    ////////////////////////////////////////////////////////////////////////////////////////////////////////
                    std::shared_ptr<QPSolver<PetscMatrix, PetscVector>> qp_solver;
                    if(use_mprgp_) {
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
                    auto box = make_lower_bound_constraints(make_ref(lb_));
                    solver.set_box_constraints(box);
                    in.get("solver", solver);
                    solver.solve(fe_problem, solution_);


                    // auto qp_solver = std::make_shared<utopia::SteihaugToint<PetscMatrix, PetscVector> >();
                    // // auto qp_solver = std::make_shared<utopia::Lanczos<PetscMatrix, PetscVector> >();
                    // qp_solver->pc_type("bjacobi"); 
                    // qp_solver->max_it(1000); 


                    // TrustRegion<PetscMatrix, PetscVector> solver(qp_solver);
                    // solver.verbose(true); 
                    // solver.delta0(1e4); 
                    // in.get("solver", solver);
                    // solver.solve(fe_problem, solution_);


                    ////////////////////////////////////////////////////////////////////////////////////////////////////////


                    auto sol_status = solver.solution_status(); 
                    const auto conv_reason = sol_status.reason;                     

                    if(adjust_dt_on_failure_ && conv_reason < 0){
                        // reset solution
                        fe_problem.get_old_solution(solution_); 

                        time_ -= dt_; 
                        dt_ = dt_ * shrinking_factor_; 
                        time_ += dt_; 
                    }
                    else{   
                        rename("X", solution_);

                        if(pressure0_!= 0.0){
                            space_.write(output_path_+"_"+std::to_string(1e-5*time_)+".vtk", solution_);      
                        }   
                        else{
                            space_.write(output_path_+"_"+std::to_string(time_)+".vtk", solution_);      
                        }

                        fe_problem.old_solution(solution_); 

                        // increment time step 
                        time_+=dt_;
                    }

                }                

            }


        private:
            FunctionSpace & space_; 
            InitialCondition<FunctionSpace> & IC_; 
            BCSetup<FunctionSpace> & BC_; 
            Scalar dt_, dt0_; 
            Scalar num_time_steps_; 
            Scalar time_; 
            std::string output_path_; 
            bool use_mprgp_; 
            bool adjust_dt_on_failure_; 
            Scalar shrinking_factor_; 
            Scalar pressure0_; 
            Scalar pressure_increase_factor_; 

            Vector solution_;
            Vector lb_; // this is quite particular for PF-frac 
    }; 



}


#endif
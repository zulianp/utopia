#ifndef UTOPIA_DM_RMTR_SETUP_HPP
#define UTOPIA_DM_RMTR_SETUP_HPP

#include "utopia_Input.hpp"
#include "utopia_Multigrid.hpp"
#include "utopia_IPTransfer.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_Multilevel.hpp"
#include "utopia_MLIncrementalLoading.hpp"
#include "utopia_PFMassMatrix.hpp"

#include <memory>

namespace utopia {

    //FIXME complete the overriding process
    template<class FunctionSpace, class ProblemType, class BCType, class ICType>
    class MLIncrementalLoading final : public IncrementalLoadingBase<FunctionSpace> {
    public:
        using Matrix   = typename FunctionSpace::Matrix;
        using Vector   = typename FunctionSpace::Vector;
        using Scalar   = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;


        MLIncrementalLoading(FunctionSpace &space_coarse, const SizeType & n_levels) : n_levels_(n_levels){
            init_ml_setup(space_coarse); 
        }   

        void read(Input &in) override {

            IncrementalLoadingBase<FunctionSpace>::read(in); 

            for (auto l=0; l < level_functions_.size(); l++){
                level_functions_[l]->read(in);    
                BC_conditions_[l]->read(in); 
            }

            IC_->read(in); 

        }


        bool init_ml_setup(FunctionSpace &space){
            return init_ml_setup(make_ref(space));
        }

        bool init_ml_setup(const std::shared_ptr<FunctionSpace> &space)
        {
            if(n_levels_ < 2) {
                std::cerr << "n_levels must be at least 2" << std::endl;
                return false;
            }

            spaces_.resize(n_levels_);
            spaces_[0] = space;

            level_functions_.resize(n_levels_); 
            auto fun = std::make_shared<ProblemType>(*spaces_[0]);
            fun->use_crack_set_irreversibiblity(false); 
            level_functions_[0] = fun; 



            BC_conditions_.resize(n_levels_); 
            auto bc = std::make_shared<BCType>(*spaces_[0]);
            BC_conditions_[0] = bc; 

            transfers_.resize(n_levels_ - 1);

            for(SizeType i = 1; i < n_levels_; ++i) {
                spaces_[i] = spaces_[i-1]->uniform_refine();

                auto fun = std::make_shared<ProblemType>(*spaces_[i]);
                

                if(i <n_levels_-1){
                    fun->use_crack_set_irreversibiblity(false); 
                }
                else{
                    fun->use_crack_set_irreversibiblity(true); 
                }


                level_functions_[i] = fun; 

                auto bc = std::make_shared<BCType>(*spaces_[i]);
                BC_conditions_[i] = bc; 


                auto I = std::make_shared<Matrix>();
                spaces_[i-1]->create_interpolation(*spaces_[i], *I);
                assert(!empty(*I));

                Matrix Iu; // = *I; 
                // MatConvert(raw_type(*I),  MATMPIAIJ, MAT_INITIAL_MATRIX, &raw_type(Iu));
                MatConvert(raw_type(*I),  I->type_override(), MAT_INITIAL_MATRIX, &raw_type(Iu));
                Matrix R = transpose(Iu); 


                PFMassMatrix<FunctionSpace> mass_matrix_assembler_fine(*spaces_[i]); 
                Matrix M_fine; 
                mass_matrix_assembler_fine.mass_matrix(M_fine); 
                // disp(H); 


                PFMassMatrix<FunctionSpace> mass_matrix_assembler_coarse(*spaces_[i-1]); 
                Matrix M_coarse; 
                mass_matrix_assembler_coarse.mass_matrix(M_coarse);                 

                Matrix inv_lumped_mass = diag(1./sum(M_coarse, 1)); 
                // disp(inv_lumped_mass); 

                Matrix P = inv_lumped_mass * R * M_fine; 
                // disp(P);
                // exit(0);

                

                // TODO:: assemble P correctly => not I^T.... 
                // transfers_[i-1] = std::make_shared<IPTransfer<Matrix, Vector> >(std::make_shared<Matrix>(Iu));    // still seq. faults 
                transfers_[i-1] = std::make_shared<MatrixTransfer<Matrix, Vector> >( std::make_shared<Matrix>(Iu), std::make_shared<Matrix>(R), std::make_shared<Matrix>(P));
            }

            


            // only needed for the finest level 
            SizeType pf_comp = 0; 
            IC_ = std::make_shared<ICType>(*spaces_.back(), pf_comp);

            return true;
        }

        FunctionSpace &fine_space()
        {
            return *spaces_.back();
        }

        const FunctionSpace &fine_space() const
        {
            return *spaces_.back();
        }

        std::shared_ptr<FunctionSpace> fine_space_ptr()
        {
            return spaces_.back();
        }

        std::shared_ptr<const FunctionSpace> fine_space_ptr() const
        {
            return spaces_.back();
        }


        ////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////

        void init_solution() override{

            spaces_.back()->create_vector(this->solution_);
            spaces_.back()->create_vector(this->lb_);
            rename("X", this->solution_);

            IC_->init(this->solution_); 
            // if(this->use_pressure_ && !this->use_constant_pressure_){
            //     Vector & pressure_vec =  level_functions_.back()->pressure_field(); 
            //     // IC_->init(this->solution_, pressure_vec); 
            // }        

            spaces_.back()->apply_constraints(this->solution_);           


            if(ProblemType * fun_finest = dynamic_cast<ProblemType *>(level_functions_.back().get())){         
                fun_finest->old_solution(this->solution_); 
            }

            // adding sol to all levels 
            for(auto l=n_levels_-1; l > 0; l--){

                ProblemType * fun_fine = dynamic_cast<ProblemType *>(level_functions_[l].get());                             
                Vector & fine_sol  = fun_fine->old_solution(); 

                ProblemType * fun_coarse = dynamic_cast<ProblemType *>(level_functions_[l-1].get());                             
                Vector & coarse_sol  = fun_coarse->old_solution();      
                spaces_[l]->create_vector(coarse_sol); 

                transfers_[l-1]->project_down(fine_sol, coarse_sol); 
                spaces_[l]->apply_constraints(coarse_sol);   

                // transfers_[l]->restrict(fine_sol, coarse_sol); 
            }

        }

        void write_to_file(FunctionSpace & space, const Scalar & time) override
        {
            // only finest level 
            IncrementalLoadingBase<FunctionSpace>::write_to_file(space, time); 

            // all levels
            // for(auto l=0; l < spaces_.size(); l++){

            //     ProblemType * fun = dynamic_cast<ProblemType *>(level_functions_[l].get());                             
            //     Vector & sol  = fun->old_solution();   
            //     rename("X", sol);

            //     spaces_[l]->write(this->output_path_+"_l_"+ std::to_string(l)+"_"+std::to_string(time)+".vtk", sol);     
            // }        

        }


        void prepare_for_solve() override{

            for(auto l=0; l < BC_conditions_.size(); l++){
                BC_conditions_[l]->emplace_time_dependent_BC(this->time_); 
            }

            // update fine level solution  and constraint 
            spaces_.back()->apply_constraints(this->solution_);    

            if(ProblemType * fun_finest = dynamic_cast<ProblemType *>(level_functions_.back().get())){    
                fun_finest->build_irreversility_constraint(this->lb_); 
            }


            for(auto l=0; l < BC_conditions_.size(); l++){
                Vector & bc_flgs    = level_functions_[l]->get_eq_constrains_flg(); 
                Vector & bc_values  = level_functions_[l]->get_eq_constrains_values(); 

                spaces_[l]->apply_constraints(bc_values);    
                spaces_[l]->build_constraints_markers(bc_flgs); 


                if(l == (BC_conditions_.size()-1)){
                    if(ProblemType * fun_finest = dynamic_cast<ProblemType *>(level_functions_.back().get())){    
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
                    std::cout<<"----- yes, constant pressure, "<< press_ts << " ......... \n";  

                    for(auto l=0; l < n_levels_; l++){
                        ProblemType * fun = dynamic_cast<ProblemType *>(level_functions_[l].get());                             
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

        void update_time_step(const SizeType & conv_reason) override
        {
            if(this->adjust_dt_on_failure_ && conv_reason < 0){
                    
                    if(ProblemType * fun_finest = dynamic_cast<ProblemType *>(level_functions_.back().get())){    
                        fun_finest->get_old_solution(this->solution_); 
                    }

                    // reset sol on all levels - important for BC conditions mostly s
                    for(auto l=n_levels_-1; l > 0; l--){

                        ProblemType * fun_fine = dynamic_cast<ProblemType *>(level_functions_[l].get());                             
                        Vector & fine_sol  = fun_fine->old_solution(); 

                        ProblemType * fun_coarse = dynamic_cast<ProblemType *>(level_functions_[l-1].get());                             
                        Vector & coarse_sol  = fun_coarse->old_solution();      
                        spaces_[l]->create_vector(coarse_sol); 

                        transfers_[l-1]->project_down(fine_sol, coarse_sol); 
                        spaces_[l]->apply_constraints(coarse_sol); 
                        // transfers_[l]->restrict(fine_sol, coarse_sol); 
                    }                    

                    this->time_ -= this->dt_; 
                    this->dt_ = this->dt_ * this->shrinking_factor_; 
                    this->time_ += this->dt_; 
                }
                else{   
                    
                    // std::cout<<"------- yes, updating...  \n"; 

                    if(ProblemType * fun_finest = dynamic_cast<ProblemType *>(level_functions_.back().get())){    
                        fun_finest->set_old_solution(this->solution_); 
                    }    

                    // update sol on all levels 
                    for(auto l=n_levels_-1; l > 0; l--){

                        ProblemType * fun_fine = dynamic_cast<ProblemType *>(level_functions_[l].get());                             
                        Vector & fine_sol  = fun_fine->old_solution(); 

                        ProblemType * fun_coarse = dynamic_cast<ProblemType *>(level_functions_[l-1].get());                             
                        Vector & coarse_sol  = fun_coarse->old_solution();      
                        spaces_[l]->create_vector(coarse_sol); 

                        transfers_[l-1]->project_down(fine_sol, coarse_sol); 
                        spaces_[l]->apply_constraints(coarse_sol); 
                        // transfers_[l]->restrict(fine_sol, coarse_sol); 
                    }     

                    if(this->pressure0_!= 0.0){
                        this->write_to_file(*spaces_.back(), 1e-5*this->time_);  
                    }   
                    else{
                        this->write_to_file(*spaces_.back(), this->time_);   
                    }

                    // increment time step 
                    this->time_ += this->dt_;
                }
        }



        void run(Input &in) override{   

            // init fine level spaces 
            this->init(in, *spaces_[n_levels_ - 1]); 



            for (auto t=1; t < this->num_time_steps_; t++)
            {
                if(mpi_world_rank()==0){
                    std::cout<<"###################################################################### \n"; 
                    std::cout<<"Time-step: "<< t << "  time:  "<< this->time_ << "  dt:  "<< this->dt_ << " \n"; 
                    std::cout<<"###################################################################### \n"; 
                }
         
                
                prepare_for_solve(); 


                // ////////////////////////////////////////////////////////////////////////////////////////////////////////
                // auto tr_strategy_fine = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
                auto tr_strategy_fine = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();
                tr_strategy_fine->atol(1e-10);
                tr_strategy_fine->verbose(false);


                auto tr_strategy_coarse = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
                // auto tr_strategy_coarse = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();
                tr_strategy_coarse->atol(1e-10);
                tr_strategy_fine->max_it(1000); 

                // TODO:: test different types of constraints
                // auto rmtr = std::make_shared<RMTR_inf<Matrix, Vector, TRGrattonBoxKornhuber<Matrix, Vector>, SECOND_ORDER> >(n_levels_);
                auto rmtr = std::make_shared<RMTR_inf<Matrix, Vector, TRKornhuberBoxKornhuber<Matrix, Vector>, GALERKIN> >(n_levels_);

                // Set TR-QP strategies
                rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
                // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

                rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
                rmtr->set_fine_tr_strategy(tr_strategy_fine);


                rmtr->set_transfer_operators(transfers_);
                rmtr->set_functions(level_functions_);                
                // rmtr->delta0(1);


                // auto box = make_lower_bound_constraints(make_ref(this->lb_));
                // rmtr->set_box_constraints(box);
                rmtr->verbose(true); 
                in.get("solver", *rmtr);


                rmtr->solve(this->solution_); 
                auto sol_status = rmtr->solution_status(); 

                //////////////////////////////////////////////////////////////////////////////////////////////////////////
                const auto conv_reason = sol_status.reason;                     
                update_time_step(conv_reason); 

            } 




        }






    private:
        SizeType n_levels_; 

        std::vector<std::shared_ptr<FunctionSpace>> spaces_;
        std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > transfers_;

        std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >  level_functions_;
        std::vector<std::shared_ptr<BCType > >  BC_conditions_;

        std::shared_ptr<ICType > IC_; 


    };

}

#endif //UTOPIA_DM_RMTR_SETUP_HPP

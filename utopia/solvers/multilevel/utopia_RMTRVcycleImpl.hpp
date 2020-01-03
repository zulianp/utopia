#ifndef UTOPIA_RMTR_BASE_IMPL_HPP
#define UTOPIA_RMTR_BASE_IMPL_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_LevelMemory.hpp"

namespace utopia
{
    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL>
    bool RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>::solve(Vector &x_h)
    {
        if(!this->check_initialization()){
            return false;
        }

        bool converged = false;
        SizeType fine_level = this->n_levels()-1;

        //-------------- INITIALIZATIONS ---------------
        this->init_memory();
        this->memory_.x[fine_level] = x_h;

        if(!this->skip_BC_checks()){
            this->make_iterate_feasible(this->function(fine_level), this->memory_.x[fine_level]);
            this->handle_equality_constraints();    
        }

        this->memory_.energy[fine_level] = this->get_multilevel_gradient_energy(this->function(fine_level),  this->memory_.s[fine_level], fine_level); 
        this->memory_.gnorm[fine_level] = this->criticality_measure(fine_level);
        this->_it_global = 0;

        //----------------------------------------------
        if(this->verbosity_level() >= VERBOSITY_LEVEL_NORMAL && mpi_world_rank() == 0)
        {
            std::cout << this->red_;
            std::string name_id = this->name() + "     Number of levels: " + std::to_string(fine_level+1)  + "   \n Fine level local dofs: " + std::to_string(this->local_level_dofs_.back());
            this->init_solver(name_id, {" it. ", "|| g ||", "   E "});

            PrintInfo::print_iter_status(this->_it_global, {this->memory_.gnorm[fine_level], this->memory_.energy[fine_level]});
            std::cout << this->def_;
        }


        while(!converged)
        {
            if(this->cycle_type() == MULTIPLICATIVE_CYCLE)
                this->multiplicative_cycle(fine_level);
            else{
                std::cout<<"ERROR::UTOPIA_RMTR << unknown cycle type, solving in multiplicative manner ... \n";
                this->multiplicative_cycle(fine_level);
            }

            #ifdef CHECK_NUM_PRECISION_mode
                if(has_nan_or_inf(this->memory_.x[fine_level]) == true)
                {
                    this->memory_.x[fine_level].set(0.0); 
                    return false;
                }
            #endif

            this->memory_.energy[fine_level] = this->get_multilevel_gradient_energy(this->function(fine_level),  this->memory_.s[fine_level], fine_level);                     
            this->memory_.gnorm[fine_level] = this->criticality_measure(fine_level);
            this->_it_global++;

            if(this->verbose() && mpi_world_rank() == 0)
            {
                std::cout << this->red_;
                if(this->verbosity_level() > VERBOSITY_LEVEL_NORMAL){
                    this->print_init_message("RMTR OUTER SOLVE", {" it. ", "|| g ||",  "   E "});
                }

                PrintInfo::print_iter_status(this->_it_global, {this->memory_.gnorm[fine_level],  this->memory_.energy[fine_level]});
                std::cout << this->def_;
            }

            // check convergence
            converged = this->check_global_convergence(this->_it_global, this->memory_.gnorm[fine_level], 9e9, this->memory_.delta[fine_level]);
        }

        // benchmarking
        NonlinearMultiLevelBase<Matrix, Vector>::print_statistics(this->_it_global);

        // passing result outside
        x_h = this->memory_.x[fine_level];
        return true;
    }











}

#endif //UTOPIA_RMTR_BASE_IMPL_HPP
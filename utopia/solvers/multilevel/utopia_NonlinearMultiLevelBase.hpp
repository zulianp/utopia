/*
* @Author: alenakopanicakova
* @Date:   2016-04-17
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-06-14
*/

#ifndef UTOPIA_NONLINEAR_ML_BASE_HPP
#define UTOPIA_NONLINEAR_ML_BASE_HPP
#include "utopia_Level.hpp"
#include "utopia_Transfer.hpp"
#include "utopia_MultiLevelBase.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include <algorithm>
#include <vector>


  namespace utopia 
  {
    /**
     * @brief      Base class for all nonlinear multilevel solvers. \n
     *             Takes care of inializing multilevel hierarchy - calls into assembly routines on each level. \n
     *             Different levels are created by interpolation and restriction.\n
     *
     * @tparam     Matrix 
     * @tparam     Vector 
     */
    template<class Matrix, class Vector, class FunctionType>
    class NonlinearMultiLevelBase : public MultiLevelBase<Matrix, Vector>,
                                    public Monitor<Matrix, Vector>
    {
    
    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Transfer<Matrix, Vector> Transfer;

      NonlinearMultiLevelBase(const Parameters params = Parameters())
      {
        set_parameters(params); 
      }

      virtual ~NonlinearMultiLevelBase(){}


      /**
       * @brief      Sets the parameters.
       *
       * @param[in]  params  The parameters
       */
      virtual void set_parameters(const Parameters params) override
      {
        MultiLevelBase<Matrix, Vector>::set_parameters(params); 

        atol_               = params.atol();            
        rtol_               = params.rtol(); 
        stol_               = params.stol(); 

        max_it_             = params.max_it(); 
        verbose_            = params.verbose(); 
        time_statistics_    = params.time_statistics();  
      }


        /**
         * @brief      The solve function for nonlinear multilevel solvers. 
         *
         * @param[in]  rhs   Function to be minimized.
         * @param[in]  rhs   The right hand side.
         * @param      x_h   The initial guess. 
         *
         */
        virtual bool solve(FunctionType &fine_fun, Vector & x_h, const Vector & rhs)
        {
            this->init_solver(this->name_id(), {" it. ", "|| grad ||", "r_norm" , "Energy"}); 
  
            bool converged = false; 
            SizeType it = 0, l = this->num_levels(); 
            Scalar r_norm, r0_norm=1, rel_norm=1, energy;

            std::cout<<"Number of levels: "<< l << "  \n"; 

            Vector g  = local_zeros(local_size(x_h)); 
            fine_fun.gradient(x_h, g); 
            r0_norm = norm2(g); 

            fine_fun.value(x_h, energy); 

            if(this->verbose())
              PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy}); 
            
            it++; 

            while(!converged)
            {            
                if(this->cycle_type() =="multiplicative")
                    this->multiplicative_cycle(fine_fun, x_h, rhs, l); 
                else if(this->cycle_type() =="full")
                {
                  this->full_cycle(fine_fun, x_h, rhs, l); 
                  this->cycle_type("multiplicative"); 
                }
                else
                    std::cout<<"ERROR::UTOPIA_Multilevel<< unknown MG type... \n"; 


                #ifdef CHECK_NUM_PRECISION_mode
                    if(has_nan_or_inf(x_h) == 1)
                    {
                        x_h = local_zeros(local_size(x_h));
                        return true; 
                    }
                #endif    

                fine_fun.gradient(x_h, g); 
                fine_fun.value(x_h, energy); 
                
                r_norm = norm2(g);
                rel_norm = r_norm/r0_norm; 

                // print iteration status on every iteration 
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy}); 

                // check convergence and print interation info
                converged = this->check_convergence(it, r_norm, rel_norm, 1); 
                it++; 
            }
            return true; 
        }

      /**
       * @brief      The solve function for nonlinear multilevel solvers. 
       *
       * @param[in]  rhs   Function to be minimized.
       * @param      x_h   The initial guess. 
       *
       */
      virtual bool solve(FunctionType & fine_fun, Vector &x_h)
      {
        Vector rhs = local_zeros(local_size(x_h)); 
        return this->solve(fine_fun,  x_h, rhs); 
      }


      /**
       * @brief      Fnction inits functions associated with assemble on each level. 
       *
       * @param[in]  level_functions  The level functions
       * @param      type             The type of ordering of level functions
       *
       */
      virtual bool init_level_functions(std::vector<FunctionType> level_functions, std::string const &type = "coarse_to_fine")
      {
          _nonlinear_levels.clear();

          if(!type.compare("fine_to_coarse"))
          {
            for(auto I = level_functions.rbegin(); I != level_functions.rend() ; ++I )
              _nonlinear_levels.push_back(std::move(*I));
          }
          else
          {
            for(auto I = level_functions.begin(); I != level_functions.end() ; ++I )
              _nonlinear_levels.push_back(std::move(*I));
          }
          return true; 
      }



       /* @brief 
                Function initializes projections  operators. 
                Operators need to be ordered FROM COARSE TO FINE. 
                If u have them in reverse order use "fine_to_coarse" flg 
                
       *
       * @param[in]  operators                The restriction operators.
       * @param      type                     Ordering of the comming operators. 
       *
       */
      virtual bool init_nonlinear_transfer(std::vector<Matrix> restriction_operators, std::vector<Matrix> projection_operators,  std::string const &type)
      {
          this->_num_levels = restriction_operators.size() + 1; 
          this->_transfers.clear();

          if(!type.compare("fine_to_coarse"))
          {
            for(auto I = restriction_operators.rbegin(), P = projection_operators.rbegin(); I != restriction_operators.rend() && P != projection_operators.rend(); ++I, ++P )
              this->_transfers.push_back(std::move(Transfer(*I, *P)));
          }
          else
          {
            for(auto I = restriction_operators.begin(), P = projection_operators.begin(); I != restriction_operators.end() &&  P != projection_operators.end() ; ++I, ++P )
              this->_transfers.push_back(std::move(Transfer(*I, *P)));
          }
          return true; 
      }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Scalar      atol() const               { return atol_; } 
        Scalar      rtol()  const              { return rtol_; } 
        Scalar      stol()  const              { return stol_; } 
        SizeType    max_it()  const            { return max_it_; } 
        bool        verbose() const                     { return verbose_; } 
        bool        time_statistics() const       { return time_statistics_; } 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void atol(const Scalar & atol_in ) { atol_ = atol_in; }; 
        void rtol(const Scalar & rtol_in ) { rtol_ = rtol_in; }; 
        void stol(const Scalar & stol_in ) { stol_ = stol_in; }; 
        void max_it(const SizeType & max_it_in ) { max_it_ = max_it_in; }; 
        void verbose(const bool & verbose_in ) { verbose_ = verbose_in; }; 
        void time_statistics(const bool & time_statistics_in ) { time_statistics_ = time_statistics_in; }; 



protected:
          /**
         * @brief      Initialization of nonlinear solver. Includes nice printout and starts calculating time of solve process. 
         *
         * @param[in]  method            The method.
         * @param[in]  status_variables  The status variables.
         */
        virtual void init_solver(const std::string &method, const std::vector<std::string> status_variables) override
        {
            if(mpi_world_rank() == 0 && verbose_)
                PrintInfo::print_init(method, status_variables); 
            
            _time.start();
        }     


        /**
         * @brief      Exit of solver. 
         *
         * @param[in]  num_it              The number iterator
         * @param[in]  convergence_reason  The convergence reason
         */
        virtual void exit_solver(const SizeType &num_it, const Scalar & convergence_reason) override
         {            
            _time.stop();
            params_.convergence_reason(convergence_reason);
            params_.num_it(num_it);

            if(verbose_)
            {
                ConvergenceReason::exitMessage_nonlinear(num_it, convergence_reason);
                if(mpi_world_rank() == 0)
                    std::cout<<"  Walltime of solve: " << _time.get_seconds() << " seconds. \n";
            }
         }


        virtual bool solver_monitor(const SizeType& it, Vector & x, Matrix & H) override
        {
          std::cout<<"utopia::NonlinearMultilevelBase:: WE ARE NOT SUPPORTING this function at the moment... \n"; 
          return true; 
        }

         /**
          * @brief      General function to check convergence in nonlinear solvers. It checks absolute, relative norm of gradient
          *             and lenght of the step size.   
          *
          * @param[in]  g_norm  The norm of the gradient. 
          * @param[in]  r_norm  The relative norm of the gradient. 
          * @param[in]  s_norm  The size of step. 
          * @param[in]  it      The number of iterations. 
          */
        virtual bool check_convergence(const SizeType &it, const Scalar & g_norm, const Scalar & r_norm, const Scalar & s_norm) override
        {   

            // termination because norm of grad is down
            if(g_norm < atol_)
            {
              exit_solver(it, ConvergenceReason::CONVERGED_FNORM_ABS);
              return true; 
            }

            // step size so small that we rather exit than wait for nan's
            if(s_norm < stol_)
            {
              exit_solver(it, ConvergenceReason::CONVERGED_SNORM_RELATIVE);
              return true; 
            }

            // step size so small that we rather exit than wait for nan's
            if(r_norm < rtol_)
            {
              exit_solver(it, ConvergenceReason::CONVERGED_FNORM_RELATIVE);
              return true; 
            }

            // check number of iterations
            if( it > max_it_)
            {
              exit_solver(it, ConvergenceReason::DIVERGED_MAX_IT);
              return true; 
            }

            return false; 
        }


        /**
         * @brief      Function looks up for ids, where we should apply Dirichlet BC and set value to required one
         *
         * @param      fun   The fun
         * @param      x     
         *
         */
        virtual bool make_iterate_feasible(FunctionType & fun, Vector & x)
        {

          Vector bc_values; 
          fun.get_boundary_values(bc_values); 

          Vector bc_ids; 
          fun.get_boundary_ids(bc_ids); 

            {
                Write<Vector> w(x);
                Read<Vector>  r_id(bc_ids);
                Read<Vector>  r_val(bc_values);

                Range range_w = range(x);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) 
                {
                    Scalar id = bc_ids.get(i);
                    Scalar value = bc_values.get(i);
                    
                    if(id == 1)
                    {
                      x.set(i, value);
                    }
                }
            }
          
          return true; 
        }


        /**
         * @brief      Function zeors correction, where we have Dirichlet BC aplied.
         *
         * @param      fun   The fun
         * @param      c     The correction
         */
        virtual bool zero_boundary_correction(FunctionType & fun, Vector & c)
        {
          Vector bc; 
          fun.get_boundary_ids(bc); 

          {
            Write<Vector> w(c);
            Read<Vector> r(bc);

            Range range_w = range(c);
            for (SizeType i = range_w.begin(); i != range_w.end(); i++) 
            {
                Scalar value = bc.get(i);
                if(value == 1)
                  c.set(i, 0);
            }
          }
          return true; 
        }



        // TODO:: make nicer => does not need to be done every iteration
        virtual bool zero_boundary_correction_mat(FunctionType & fun, Matrix & M)
        {
          Vector bc; 
          fun.get_boundary_ids(bc); 
          std::vector<SizeType> index;

          {
            Read<Vector> r(bc);

            Range range_w = range(bc);
            for (SizeType i = range_w.begin(); i != range_w.end(); i++) 
            {
                if(bc.get(i) == 1)
                  index.push_back(i); 
            }
          }
          set_zero_rows(M, index); 
          return true; 
        }


        inline FunctionType &levels(const SizeType &l) { return this->_nonlinear_levels[l];  }
        inline Transfer &transfers(const SizeType & l) { return this->_transfers[l];  }


        /**
         * @brief     Multiplicative V/W cycle
         *
         * @param      fine_fun  Function to be minimized 
         * @param      u_l       The iterate 
         * @param[in]  f         Right hand side 
         * @param[in]  l         level
         *                       
         */
        virtual bool multiplicative_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & l)=0; 

        /**
         * @return     Name of solver - to have nice printouts 
         */
        virtual std::string name_id() = 0; 


        /**
         * @brief      Coarse solve function. 
         *
         * @param      fun   Function to be minimized 
         * @param      x     Iterate
         * @param[in]  rhs   The right hand side
         *
         */
        virtual bool coarse_solve(FunctionType &fun, Vector &x, const Vector & rhs) = 0; 


        /**
         * @brief     Full multigrid cycle, after running F cycle once and reaching discretization tolerance, it continues as multiplicative cycle to achieve prescribed tolerance
         *
         * @param      fine_fun  Function to be minimized 
         * @param      u_l       The iterate 
         * @param[in]  f         Right hand side 
         * @param[in]  l         level
         *                       
         */
        virtual bool full_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & l)
        {            
          for(SizeType i = l-2; i >=0; i--)
          {
            transfers(i).restrict(u_l, u_l); 
            this->make_iterate_feasible(levels(i), u_l); 
          }

          Vector L_l = local_zeros(local_size(u_l));
          this->coarse_solve(levels(0), u_l, L_l); 

          transfers(0).interpolate(u_l, u_l); 

          for(SizeType i = 1; i <l-1; i++)
          {
            for(SizeType j = 0; j < this->v_cycle_repetition(); j++)
            {
              Vector f = local_zeros(local_size(u_l));
              this->multiplicative_cycle(levels(i), u_l, f, i+1); 
            }
            transfers(i).interpolate(u_l, u_l); 
          }
            return true; 
        }



    protected:
        std::vector<FunctionType>                      _nonlinear_levels;  
        Parameters params_;                           /*!< Solver parameters. */  

        
        // ... GENERAL SOLVER PARAMETERS ...
        Scalar atol_;                   /*!< Absolute tolerance. */  
        Scalar rtol_;                   /*!< Relative tolerance. */  
        Scalar stol_;                   /*!< Step tolerance. */  

        SizeType max_it_;               /*!< Maximum number of iterations. */  
        SizeType verbose_;              /*!< Verobse enable? . */  
        SizeType time_statistics_;      /*!< Perform time stats or not? */  


        Chrono _time;                 /*!<Timing of solver. */

  };

}

#endif //UTOPIA_NONLINEAR_ML_BASE_HPP


/*
* @Author: alenakopanicakova
* @Date:   2016-04-17
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-05-10
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
     *             Takes care of inializing multilevel hierarchy. \n
     *             Different levels are created by interpolation and restriction operators.\n
     *             Additionally, it provides stifness matrices on each level, created by using Galerkin assembly. \n
     *
     * @tparam     Matrix 
     * @tparam     Vector 
     */
    template<class Matrix, class Vector, class FunctionType>
    class NonlinearMultiLevelBase : public MultiLevelBase<Matrix, Vector>,
                                    public Monitor<Matrix, Vector>
    {

      // typedef utopia::Function<Matrix, Vector> Function;
    
    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Transfer<Matrix, Vector> Transfer;

      NonlinearMultiLevelBase(const Parameters params = Parameters())
      {
        set_parameters(params); 
      }

      virtual ~NonlinearMultiLevelBase(){}


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
          std::cout<<"utopia::NonlinearMultilevelBase:: WE ARE NOT SUPPORTING monitor for FAS at the moment... \n"; 
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





        virtual bool make_iterate_feasible(FunctionType & fun, Vector & x)
        {

          // std::cout<<"make_iterate_feasible   \n"; 
          Vector bc; 
          fun.get_boundary_values(bc); 

          // std::cout<<"yes non zero:   "; 

            {
                Write<Vector> w(x);
                Read<Vector> r(bc);

                Range range_w = range(x);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) 
                {
                    Scalar value = bc.get(i);
                    
                    if(value != 0)
                    {
                      // std::cout<<"   "<< i << "    "; 
                      x.set(i, value);
                    }
                }
            }

            std::cout<<"    \n"; 
          return true; 
        }



        // zero correction
        virtual bool zero_correction_contributions(FunctionType & fun, Vector & c)
        {

          // std::cout<<"zero_correction_contributions   \n"; 
          Vector bc; 
          fun.get_boundary_values(bc); 


            {
                Write<Vector> w(c);
                Read<Vector> r(bc);

                Range range_w = range(c);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) 
                {
                    Scalar value = bc.get(i);
                    
                    if(value != 0)
                    {
                      // std::cout<<"yes non zero: "<< i << "   \n"; 
                      c.set(i, 0);
                    }
                }
            }

          

          return true; 
        }



    protected:
        std::vector<FunctionType>                      _nonlinear_levels;  
        Parameters params_;                /*!< Solver parameters. */  



       // std::vector<Transfer>               _transfers;   /*!< vector of transfer operators  */
        
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


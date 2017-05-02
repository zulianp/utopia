/*
* @Author: alenakopanicakova
* @Date:   2016-03-28
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-05-01
*/

#ifndef UTOPIA_ML_BASE_HPP
#define UTOPIA_ML_BASE_HPP
#include "utopia_Level.hpp"
#include "utopia_Transfer.hpp"
#include "utopia_Core.hpp"
#include <algorithm>
#include <vector>



  namespace utopia 
  {
    /**
     * @brief      Base class for all multilevel solvers. \n
     *             Takes care of inializing multilevel hierarchy. \n
     *             Different levels are created by interpolation and restriction operators.\n
     *             Additionally, it provides stifness matrices on each level, created by using Galerkin assembly. \n
     *
     * @tparam     Matrix 
     * @tparam     Vector 
     */
    template<class Matrix, class Vector>
    class MultiLevelBase
    {
      typedef UTOPIA_SCALAR(Vector)    Scalar;
      typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
      typedef utopia::Level<Matrix, Vector> Level;
      typedef utopia::Transfer<Matrix, Vector> Transfer;

    
    public:

      MultiLevelBase(const Parameters params = Parameters())
      {
        set_parameters(params); 
      }

      virtual ~MultiLevelBase(){}

      virtual void set_parameters(const Parameters params)
      {
          _parameters = params; 
      
          _pre_smoothing_steps = params.pre_smoothing_steps(); 
          _post_smoothing_steps = params.post_smoothing_steps(); 
          _mg_type = params.mg_type(); 
          _cycle_type = params.cycle_type(); 
          _v_cycle_repetition = 1;  // TODO:: create option in params for this 
      }


      /**
       * @brief 
                Function initializes restriction transfer operators. 
                Operators need to be ordered FROM COARSE TO FINE. 
                If u have them in reverse order use "fine_to_coarse" flg 
                
       *
       * @param[in]  operators                The restriction operators.
       * @param      type                     Ordering of the comming operators. 
       *
       */
      virtual bool init_transfer(std::vector<Matrix> restriction_operators, std::string const &type = "coarse_to_fine")
      {
          _num_levels = restriction_operators.size() + 1; 
          _transfers.clear();

          if(!type.compare("fine_to_coarse"))
          {
            for(auto I = restriction_operators.rbegin(); I != restriction_operators.rend() ; ++I )
              _transfers.push_back(std::move(Transfer(*I)));
          }
          else
          {
            for(auto I = restriction_operators.begin(); I != restriction_operators.end() ; ++I )
              _transfers.push_back(std::move(Transfer(*I)));
          }
          return true; 
      }


      
      /**
       * @brief 
       *        The function creates corser level operators by using Galerkin assembly. 
       *        
       *        $\f J_{i-1} = R * J_{i} * I  $\f
       *        
       *        Resulting operators are assumed to go from fines = 0 to coarse = num_levels
       *             
       * @param[in]  stifness matrix for finest level
       *
       */
      virtual bool galerkin_assembly(const Matrix & A)
      {
          _levels.clear();
          SizeType t_s = _transfers.size(); 
          if(t_s <= 0)
            std::cout<<"Provide interpolation operators first!  \n"; 

          _levels.push_back(std::move(A));  
        
          for(SizeType i = 1; i < _num_levels; i++)
          {
            // J_{i-1} = R * J_{i} * I
            Matrix J_h; 
            _transfers[t_s - i].restrict(_levels[i-1].A(), J_h); 
            _levels.push_back(std::move(J_h));              
          }
          
          std::reverse(std::begin(_levels), std::end(_levels));
          return true; 
      }




      /**
       * @brief 
       *        The function creates corser level operators provided by assembling on differnet levels of MG hierarchy
       *        
       * @param[in]  stifness matrix for finest level
       *
       */
      virtual bool assembly_linear_operators(const std::vector<Matrix> A, std::string const &type = "coarse_to_fine")
      {
          _levels.clear();
          if(!type.compare("fine_to_coarse"))
          {
            for(auto I = A.rbegin(); I != A.rend() ; ++I )
              _levels.push_back(std::move(*I));
          }
          else
          {
            for(auto I = A.begin(); I != A.end() ; ++I )
              _levels.push_back(std::move(*I));
          }
          std::cout<<"size of levels: "<< _levels.size() << "  \n"; 
          return true; 
      }




        /**
         * @brief      Returns number of levels in hierarchy.
         */
        virtual SizeType num_levels()
        {
          return _num_levels; 
        }



        /**
         * @brief      Function sets type of cycle
         */
        bool cycle_type(const std::string & type_in)
        {
            _cycle_type = type_in; 
            return true; 
        }

        /**
         * @brief    Sets amount of V-cycles inside of F-cycle
         */
        bool v_cycle_repetition(const SizeType & v_cycle_repetition_in)
        {
            _v_cycle_repetition = v_cycle_repetition_in; 
            return true; 
        }


        /**
         * @brief      Setting number pre-smoothing steps. 
         *
         * @param[in]  pre_smoothing_steps_in  Number of pre-smoothing steps.
         */
        void pre_smoothing_steps(const SizeType & pre_smoothing_steps_in ) 
        { 
            _pre_smoothing_steps = pre_smoothing_steps_in; 
        }; 

        /**
         * @brief      Setting number of post-smoothing steps.
         *
         * @param[in]  post_smoothing_steps_in  Number of post-smoothing steps.
         */
        void post_smoothing_steps(const SizeType & post_smoothing_steps_in )
        { 
            _post_smoothing_steps = post_smoothing_steps_in; 
        }; 

        /**
         * @brief      Setting type of MG:  1 goes for V_CYCLE, 2 for W-cycle. 
         *
         * @param[in]  mg_type_in  Choice of MG cycle.
         */
        void mg_type(const bool & mg_type_in ) 
        { 
            _mg_type = mg_type_in; 
        }; 

        /**
         * @return     Number of pre-smoothing steps. 
         */
        SizeType  pre_smoothing_steps() const         
        { 
            return _pre_smoothing_steps; 
        } 

        /**
         * @return     Number of post-smoothing steps.
         */
        SizeType  post_smoothing_steps() const        
        { 
            return _post_smoothing_steps; 
        } 

        /**
         * @return     Type of MG cycle. 
         */
        bool  mg_type() const                     { return _mg_type; } 

        /**
         * @return     Type of MG cycle. 
         */
        std::string  cycle_type() const         { return _cycle_type; } 


        /**
         * @brief      Amount of V-cycles on each level during full-cycle
         */
        SizeType v_cycle_repetition() const {return _v_cycle_repetition; }





    protected:
        SizeType _num_levels;                             /*!< number of levels in ML   -n   */ 
        std::vector<Level>                  _levels;      /*!< vector of level operators     */
        std::vector<Transfer>               _transfers;   /*!< vector of transfer operators  */

        Parameters                          _parameters; 

        SizeType                            _pre_smoothing_steps; 
        SizeType                            _post_smoothing_steps; 
        SizeType                            _mg_type; 

        std::string                         _cycle_type; 
        SizeType                            _v_cycle_repetition; 
  };

}

#endif //UTOPIA_ML_BASE_HPP


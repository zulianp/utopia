/*
* @Author: alenakopanicakova
* @Date:   2016-03-28
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-01-29
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

      MultiLevelBase( ){}

      virtual ~MultiLevelBase(){}


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
       * @brief      Returns number of levels in hierarchy.
       */
      SizeType num_levels()
      {
        return _num_levels; 
      }


    protected:
        SizeType _num_levels;                                 /*!< number of levels in ML   -n   */ 
        std::vector<Level>                      _levels;      /*!< vector of level operators     */
        std::vector<Transfer>                   _transfers;   /*!< vector of transfer operators  */
        
        SizeType mpi_size = mpi_world_size();
        SizeType mpi_rank = mpi_world_rank();
  };

}

#endif //UTOPIA_ML_BASE_HPP


#ifndef UTOPIA_SOLVER_TRUSTREGION_BOX_BASE_HPP
#define UTOPIA_SOLVER_TRUSTREGION_BOX_BASE_HPP

#include <algorithm>

#include "utopia_TRBase.hpp"
#include "utopia_Parameters.hpp"    
#include "utopia_NumericalTollerance.hpp"
#include "utopia_BoxConstraints.hpp"

namespace utopia  
{
  template<class Matrix, class Vector>
  /**
   * @brief      Base class for  TR bound constrained solvers. Contains all general routines related to TR solvers.
   *             Design of class allows to provide different TR strategies in order to solve TR subproblem. 
   */ 
  class TrustRegionBoxBase : public TrustRegionBase<Matrix, Vector>
  {
    typedef UTOPIA_SCALAR(Vector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

    typedef utopia::BoxConstraints<Vector>              BoxConstraints;
    

  public:
      TrustRegionBoxBase( const Parameters params = Parameters())
      {
          set_parameters(params); 
      }

      virtual void set_parameters(const Parameters params) override
      {
        TrustRegionBase<Matrix, Vector>::set_parameters(params);
      }

      /**
     * @brief      Sets the box constraints.
     *
     * @param      box   The box
     *
     */
    virtual bool set_box_constraints(BoxConstraints & box)
    {
      box_constraints_ = box; 
      return true; 
    }

    /**
     * @brief      Gets the box constraints.
     *
     * @return     The box constraints.
     */
    virtual BoxConstraints & get_box_constraints()
    {
      return box_constraints_; 
    }



  protected:

      /**
       * @brief      Computes criticality measure, which is later used for termination conditions
       *
       * @param[in]  x     The current iterate
       * @param[in]  g     The gradient
       *
       */
      virtual Scalar criticality_measure_infty(const Vector & x, const Vector & g)
      {

        Vector Pc; 
        Vector x_g = x - g; 
        Vector ub, lb; 

        Scalar n = local_size(x).get(0); 

        if(box_constraints_.has_upper_bound())
          ub = *box_constraints_.upper_bound(); 
        else
          ub =  local_values(n, 9e15); 

        if(box_constraints_.has_lower_bound())
          lb = *box_constraints_.lower_bound(); 
        else
          lb =  local_values(n, -9e15); 

        get_projection(x_g, lb, ub, Pc); 
        
        Pc -= x; 
        return norm2(Pc); 
      }



      /*!
      \details
      TR radius update function
      @note
      \param rho            - ared/pred
      \param radius          - tr. radius
      \param p_k            - iterate step
      */
    virtual bool delta_update(const Scalar &rho, const Vector &p_k, Scalar &radius) override
    {
      if(rho < this->eta1())
      {
        radius = std::max(Scalar(this->gamma1() * norm_infty(p_k)), this->delta_min()); 
      }
      else if (rho > this->eta2() )
      {
        // Scalar intermediate = std::max(Scalar(this->gamma2() * norm_infty(p_k)), radius); 

        Scalar intermediate = this->gamma2() * radius; 
        radius = std::min(intermediate, this->delta_max()); 
      }      
      return true; 
    }


      /**
       * @brief      Projection onto feasible set
       *
       * @param[in]  x     The quantity to be projected 
       * @param[in]  lb    The lower bound
       * @param[in]  ub    The upper bound
       * @param      Pc    Projection
       *
       */
      bool get_projection(const Vector & x, const Vector &lb, const Vector &ub, Vector & Pc)
      {
          Pc = local_values(local_size(x).get(0), 1.0);; 
          {
              Read<Vector> r_ub(ub), r_lb(lb), r_x(x);
              Write<Vector> wv(Pc); 

              each_write(Pc, [&ub, &lb, &x](const SizeType i) -> double { 
                          Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  
                          if(li >= xi)
                            return li; 
                          else
                            return (ui <= xi) ? ui : xi; }   );
          }

          return true;
      }

      virtual BoxConstraints  merge_tr_with_pointwise_constrains(const Vector & x_k, const Scalar & radius)
      {
          Vector l_f, u_f; 

          if(box_constraints_.has_upper_bound())
          {
              Vector u =  *box_constraints_.upper_bound() - x_k; 
              u_f = local_zeros(local_size(x_k).get(0)); 
              {   
                  Read<Vector> rv(u); 
                  Write<Vector> wv(u_f); 

                  each_write(u_f, [radius, &u](const SizeType i) -> double { 
                      const auto val = u.get(i);
                      return  (val <= radius) ? val : radius; }   );
              }
          }
          else
              u_f = radius * local_values(local_size(x_k).get(0), 1.0); ; 

          if(box_constraints_.has_lower_bound())
          {
              Vector l = *(box_constraints_.lower_bound()) - x_k; 
              l_f = local_zeros(local_size(x_k).get(0)); 

              {   
                  Read<Vector> rv(l); 
                  Write<Vector> wv(l_f); 

                  each_write(l_f, [radius, &l](const SizeType i) -> double { 
                      return  (l.get(i) >= -1*radius)  ? l.get(i) : -1 * radius;  }   );
              }
          }
          else
              l_f = -1 * radius * local_values(local_size(x_k).get(0), 1.0); 

          return make_box_constaints(std::make_shared<Vector>(l_f), std::make_shared<Vector>(u_f));
      }


    protected: 
      BoxConstraints box_constraints_; 

  };

}

#endif //UTOPIA_SOLVER_TRUSTREGION_BOX_BASE_HPP


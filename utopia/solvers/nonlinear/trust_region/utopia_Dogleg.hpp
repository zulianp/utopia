/*
* @Author: alenakopanicakova
* @Date:   2016-05-11
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
*/

#ifndef UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP
#define UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP
#include "utopia_TRSubproblem.hpp"
#include "utopia_CauchyPoint.hpp"
#include "utopia_Parameters.hpp"    
#include "utopia_LinearSolverInterfaces.hpp"

// #ifdef WITH_PETSC
//     #include "utopia_PETScFactorization.hpp"
// #endif //WITH_PETSC

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// NEEEEEDS PROPER INTERFACE ///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace utopia 
{

	template<class Matrix, class Vector>
    class Dogleg : public TRSubproblem<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;

    public:

        Dogleg(const Parameters & params = Parameters()): 
                                                          TRSubproblem<Matrix, Vector>(params)
        {

        };

protected: 
        bool unpreconditioned_solve(const Matrix &B, const Vector &g, Vector &p_k) override
        {
            Vector p_N = g, p_SD = g, p_CP = g;
            Scalar pred_N, g_B_g = dot(g, B * g), pNlen; 

            // this is the worst hard-codding ever 
            auto lsolver = std::make_shared< Factorization<Matrix, Vector> >();
			#ifdef PETSC_HAVE_MUMPS
            	lsolver->set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
			#endif //PETSC_HAVE_MUMPS
			
            //
            lsolver->solve(B, -1 * g, p_N);


            // model reduction for newton point
            get_pred(g, B, p_N, pred_N); 


            // good situation: p_N reduces m_k(p)
	    	if(pred_N >= 0)
            {   
                pNlen = norm2(p_N);
                // newton point inside of TR
                if(pNlen <= this->current_radius())
                {
                    p_k = p_N;
                }
                else
                {
                    // check if CP or SD 
                    //steepest descent step 
                    p_k = (dot(g, g)/g_B_g) * g;
                    p_k *= -1.0;
                    Scalar SD_norm = norm2(p_k); 
                  //  std::cout<<"pNlen <= SD taking   \n"; 

                    if(SD_norm < this->current_radius())
                    {
                    	//std::cout<<"pNlen <= inside  \n"; 
                        // p_SD inside, p_N outside => dogleg 
                        // p_k = p_CP + \tau (p_N - P_CP), 
                		// where \tau is the largest value in [0,1]
                        // such that || p_k || <= \delta                   	
                        Vector d = p_N - p_SD; 
                        this->quad_solver(p_CP, d, this->current_radius(),  p_k);
                    }
                }

            }

            //  bad situation: p_N increases m_k(p)
            else
            {
                if(g_B_g > 5 * std::sqrt(1e-14))
                {
                    //  1D function is convex enough
                    //  find unconstrained SD minimizer for model fun
                    p_SD = -(dot(g, g)/g_B_g) * g;
                    Scalar pSDlen = norm2(p_SD);

                    if(pSDlen <= this->current_radius())
                    {
                        p_k = p_SD;
                    }
                    else
                    {
                        p_k = this->current_radius() * (p_SD* (1/pSDlen));
                    }
                }
                else
                {
                    //  1D slice of model function in gradient direction is concave:
                    //  SD minimizer lies on TR boundary, in direction of -g'
                    cauchy_point.unpreconditioned_solve(B, g, p_k); 
                }
            }

	    	return true; 

        }

    virtual bool get_pred(const Vector & g, const Matrix & B, const Vector & p_k, Scalar &pred)
    {
      Scalar l_term = dot(g, p_k);
      Scalar qp_term = dot(p_k, B * p_k);
      pred = - l_term - 0.5 * qp_term; 
      return true; 
    }


    private: 
        CauchyPoint<Matrix, Vector> cauchy_point;       /*!< Cauchy Point  */  


    };
}


#endif //UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP


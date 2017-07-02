/*
* @Author: alenakopanicakova
* @Date:   2016-05-11
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-02
*/

#ifndef UTOPIA_TR_SUBPROBLEM_STEIHAUG_TOINT_HPP
#define UTOPIA_TR_SUBPROBLEM_STEIHAUG_TOINT_HPP
#include "utopia_TRSubproblem.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Preconditioner.hpp"

namespace utopia 
{
	/**
	 * @brief      Class for Steihaug Toint conjugate gradient.
	 */
	template<class Matrix, class Vector>
    class SteihaugToint : public TRSubproblem<Matrix, Vector>
    {
		typedef UTOPIA_SCALAR(Vector) Scalar;

    public:

    	SteihaugToint(const Parameters params = Parameters()): 
    														TRSubproblem<Matrix, Vector>(params)
        {  };

        virtual ~SteihaugToint(){}

private:
        bool unpreconditioned_solve(const Matrix &B, const Vector &g, Vector &p_k) override
        {

			Vector r = -1 * g, d = r, s = local_zeros(local_size(g)), s1 = s; 
	    	Scalar alpha, g_norm, d_B_d, z, z1, tau;
	    	SizeType it = 0;

	    	p_k = local_zeros(local_size(p_k)); 
	    	g_norm = norm2(g);
	    	z = g_norm * g_norm; 


	    	// this->atol(std::min(0.5, std::sqrt(g_norm)) * g_norm); 
	    	// this->atol(1e-14); 
	    	// this->rtol(1e-14); 
	    	// this->stol(1e-14); 

	    	this->init_solver(" Utopia Steihaug-Toint CG ", {"it. ", "||r||" }); 
            bool converged = false; 

        	while(!converged)
        	{
	    		d_B_d = dot(d, B * d); 

	    		if(d_B_d <= 0)
	    		{
	    			s = p_k;
	    			tau = this->quad_solver(s, d, this->current_radius(),  p_k);
	    			return true; 
	    		}

	    		alpha = z / d_B_d; 
	    		s1 = p_k + alpha * d;

	    		if(norm2(s1) >= this->current_radius())
	    		{
	    			s = p_k; 
	    			tau = this->quad_solver(s, d, this->current_radius(),  p_k);
	    			return true; 
	    		}
	    	
	    		p_k = s1; 
	    		r -= (alpha * B * d); 

	    		z1 = dot(r,r);
	    		d = r + (z1/z) * d;

	    		z = z1; 

        		g_norm = std::sqrt(z); 
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {g_norm}); 

                converged = this->check_convergence(it, g_norm, 1, 1); 
                it++; 
	    	}

        	return true; 
        }


        bool preconditioned_solve(const Matrix & /*B*/, const Vector &/*g*/, Vector &/*p_k*/) override
        {
        	std::cout<<"SteihaugToint:: preconditioned solve not imlemented yet ... \n"; 
        	return false; 
        }


    };


}

#endif //UTOPIA_TR_SUBPROBLEM_STEIHAUG_TOINT_HPP
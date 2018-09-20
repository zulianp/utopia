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
	template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class SteihaugToint : public TRSubproblem<Matrix, Vector>
    {
		typedef UTOPIA_SCALAR(Vector) Scalar;

    public:

    	SteihaugToint(const Parameters params = Parameters()):
    				  TRSubproblem<Matrix, Vector>(params)
        {  };

        virtual ~SteihaugToint(){}

        SteihaugToint * clone() const override
        {
        	return new SteihaugToint();
        }

	protected:
        bool unpreconditioned_solve(const Matrix &B, const Vector &g, Vector &p_k) override
        {
			Vector r = -1 * g, d = r, s = local_zeros(local_size(g)), s1 = s;
	    	Scalar alpha, g_norm, d_B_d, z, z1;
	    	SizeType it = 1;

	    	p_k = local_zeros(local_size(p_k));
	    	g_norm = norm2(g);
	    	z = g_norm * g_norm;


	    	this->init_solver(" ST-CG ", {"it. ", "||r||" });
            bool converged = false;

        	while(!converged)
        	{
	    		d_B_d = dot(d, B * d);

	    		if(d_B_d <= 0)
	    		{
	    			s = p_k;
	    			this->quad_solver(s, d, this->current_radius(),  p_k);
	    			return true;
	    		}

	    		alpha = z / d_B_d;
	    		s1 = p_k + alpha * d;

	    		if(norm2(s1) >= this->current_radius())
	    		{
	    			s = p_k;
	    			this->quad_solver(s, d, this->current_radius(),  p_k);
	    			return true;
	    		}

	    		p_k = s1;
	    		r -= alpha * (B * d);

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


        bool preconditioned_solve(const Matrix &B, const Vector &g, Vector &s_k) override
        {
        	bool converged = false;
            SizeType it=0; 

			Vector v_k = local_zeros(local_size(g)); 
				   s_k = local_zeros(local_size(g)); 
			Vector g_k = g; 

			Scalar g_norm = norm2(g_k);  

			this->init_solver(" Precond-ST-CG ", {"it. ", "||g||" });
    		if(this->verbose())
                PrintInfo::print_iter_status(it, {g_norm});
            it++; 

			this->precond_->apply(g_k, v_k);

			// if preconditioner yields nans or inf, return gradient step 
			if(has_nan_or_inf(v_k))
	    	{
	    		s_k = g_k; // TODO:: check for minus sign 
	    		return false; 
	    	}   

			Vector p_k = -1.0 * v_k; 

            Scalar alpha, kappa, betta; 
            Scalar g_v_prod_old, g_v_prod_new; 

            Scalar s_norm=0.0, s_norm_new=0.0,  sMp=0.0; 
            Scalar p_norm = dot(g_k, v_k); 
            Scalar r2 = this->current_radius() * this->current_radius(); 

        	while(!converged)
        	{
        		Vector B_p_k = B* p_k; 
	    		kappa = dot(p_k,B_p_k);

	    		// identify negative curvature 
	    		if(kappa <= 0.0)
	    		{
	    			Scalar term1 = sMp*sMp + (p_norm  * (r2 - s_norm)); 
		    		Scalar tau = (std::sqrt(term1) - sMp)/p_norm; 
		    		s_k += tau * p_k; 

	    			return true;
	    		}

	    		g_v_prod_old = dot(g_k, v_k); 
	    		alpha = g_v_prod_old/kappa; 

	    		s_norm_new = s_norm + (2.0* alpha * sMp) + (alpha * alpha * p_norm); 

	    		// ||s_k||_M > \Delta => terminate
	    		if(s_norm_new >= this->current_radius())
	    		{	
	    			Scalar term1 = sMp*sMp + (p_norm  * (r2 - s_norm)); 
		    		Scalar tau = (std::sqrt(term1) - sMp)/p_norm; 
		    		s_k += tau * p_k; 

	    			return true;
	    		}

				s_k += alpha * p_k; 	    		
	    		g_k += alpha * B_p_k; 
	    		v_k = local_zeros(local_size(g_k));

	    		this->precond_->apply(g_k, v_k);

	    		// if precond yields nans or infs, we just return current iterate 
	    		if(has_nan_or_inf(v_k))
		    		return false; 

	    		g_v_prod_new = dot(g_k, v_k); 
	    		betta  = g_v_prod_new/ g_v_prod_old; 
	    		p_k = betta * p_k - v_k; 
 		

	    		// updating norms recursively 
	    		sMp = (betta * sMp) + (alpha * p_norm); 
	    		p_norm = g_v_prod_new + (betta*betta * p_norm); 
	    		s_norm = s_norm_new; 

	    		g_norm = norm2(g_k);  

	    		if(this->verbose())
                    PrintInfo::print_iter_status(it, {g_norm});
	    		
	    		if(has_nan_or_inf(g_k))
	    		{
	    			return false; 
	    		}                

                converged = this->check_convergence(it, g_norm, 1, 1);
                it++;
	    	}

        	return true;
        }


    };


}

#endif //UTOPIA_TR_SUBPROBLEM_STEIHAUG_TOINT_HPP
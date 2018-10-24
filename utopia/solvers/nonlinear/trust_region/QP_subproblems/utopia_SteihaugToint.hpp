#ifndef UTOPIA_TR_L2_SUBPROBLEM_STEIHAUG_TOINT_HPP
#define UTOPIA_TR_L2_SUBPROBLEM_STEIHAUG_TOINT_HPP
#include "utopia_TRSubproblem.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Preconditioner.hpp"

namespace utopia
{


	/**
	 * @brief      Class for Steihaug Toint conjugate gradient.
	 */
	template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class SteihaugToint : public TRSubproblem<Matrix, Vector>, public MatrixFreeLinearSolver<Vector>
    {
		typedef UTOPIA_SCALAR(Vector) Scalar;

    public:

    	using TRSubproblem<Matrix, Vector>::tr_constrained_solve; 
    	using TRSubproblem<Matrix, Vector>::set_preconditioner; 
    	

    	SteihaugToint(const Parameters params = Parameters()):
    				  TRSubproblem<Matrix, Vector>(params)
        {  };

        virtual ~SteihaugToint(){}

        virtual SteihaugToint * clone() const override
        {
        	return new SteihaugToint(*this);
        }


        virtual bool apply(const Vector &b, Vector &x) override
        {
        	init(local_size(b).get(0)); 
            if(this->precond_) 
            {
            	auto A_ptr = utopia::op(this->get_operator());
                return preconditioned_solve(*A_ptr, b, x);
            } else {
            	auto A_ptr = utopia::op(this->get_operator());
                return unpreconditioned_solve(*A_ptr, b, x);
            }
        }

        virtual bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override 
        {
            utopia_error("ProjectedGradientActiveSet missing solve implementation.... \n"); 
            return false; 
        }

	public:
        virtual bool unpreconditioned_solve(const Operator<Vector> &B, const Vector &g, Vector &corr) override
        {
        	init(local_size(g).get(0)); 
			r = -1 * g; 
			v_k = r; 
	    	Scalar alpha, g_norm, d_B_d, z, z1;
	    	SizeType it = 1;

	    	corr = local_zeros(local_size(g));
	    	g_norm = norm2(g);
	    	z = g_norm * g_norm;


	    	this->init_solver(" ST-CG ", {"it. ", "||r||" });
            bool converged = false;

        	while(!converged)
        	{
        		B.apply(v_k, B_p_k); 
	    		d_B_d = dot(v_k, B_p_k);

	    		if(d_B_d <= 0)
	    		{
	    			Vector s = corr;
	    			this->quad_solver(s, v_k, this->current_radius(),  corr);
	    			return true;
	    		}

	    		alpha = z / d_B_d;
	    		p_k = corr + alpha * v_k;

	    		if(norm2(p_k) >= this->current_radius())
	    		{
	    			Vector s = corr;
	    			this->quad_solver(s, v_k, this->current_radius(),  corr);
	    			return true;
	    		}

	    		corr = p_k;
	    		r -= alpha * B_p_k;

	    		z1 = dot(r,r);
	    		v_k = r + (z1/z) * v_k;

	    		z = z1;

        		g_norm = std::sqrt(z);

                if(this->verbose())
                    PrintInfo::print_iter_status(it, {g_norm});

                converged = this->check_convergence(it, g_norm, 1, 1);
                it++;
	    	}

        	return true;
        }


        virtual bool preconditioned_solve(const Operator<Vector> &B, const Vector &g, Vector &s_k) override
        {
        	init(local_size(g).get(0)); 
        	bool converged = false;
            SizeType it=0; 

			s_k = local_zeros(local_size(g)); 
			r = g; 

			Scalar g_norm = norm2(r);  

			this->init_solver(" Precond-ST-CG ", {"it. ", "||g||", "||s||", "||p||", "sMp" });
    		if(this->verbose())
                PrintInfo::print_iter_status(it, {g_norm});
            it++; 

			this->precond_->apply(r, v_k);

			p_k = -1.0 * v_k; 

            Scalar alpha, kappa, betta; 
            Scalar g_v_prod_old, g_v_prod_new; 

            Scalar s_norm=0.0, s_norm_new=0.0,  sMp=0.0; 
			Scalar r2 = this->current_radius() * this->current_radius(); 


            Scalar p_norm = dot(r, v_k); 

			// if preconditioner yields nans or inf, or is precond. dir is indefinite - return gradient step 
			if(!std::isfinite(p_norm) || p_norm < 0.0)
	    	{
	    		Scalar alpha_termination; 
	    		if(r2 >= g_norm)
	    			alpha_termination = 1.0;  		// grad. step is inside of tr boundary, just take it
	    		else
	    			alpha_termination = std::sqrt(r2/g_norm);  // grad. step is outside of tr boundary, project on the boundary

	    		if(std::isfinite(alpha_termination))
	    			s_k -= alpha_termination * r;  

	    		return true; 
	    	}   


        	while(!converged)
        	{
        		B.apply(p_k, B_p_k); 
	    		kappa = dot(p_k,B_p_k);

	    		// identify negative curvature 
	    		if(kappa <= 0.0)
	    		{
	    			Scalar term1 = sMp*sMp + (p_norm  * (r2 - s_norm)); 
		    		Scalar tau = (std::sqrt(term1) - sMp)/p_norm; 

		    		if(std::isfinite(tau))
		    			s_k += tau * p_k; 

	    			return true;
	    		}

	    		g_v_prod_old = dot(r, v_k); 
	    		alpha = g_v_prod_old/kappa; 

	    		s_norm_new = s_norm + (2.0* alpha * sMp) + (alpha * alpha * p_norm); 

	    		// ||s_k||_M > \Delta => terminate
	    		if(s_norm_new >= this->current_radius())
	    		{	
	    			Scalar term1 = sMp*sMp + (p_norm  * (r2 - s_norm)); 
		    		Scalar tau = (std::sqrt(term1) - sMp)/p_norm; 

		    		if(std::isfinite(tau))
		    			s_k += tau * p_k; 

	    			return true;
	    		}

	    		if(std::isfinite(alpha))
					s_k += alpha * p_k; 	    		
				else
					return false; 

	    		r += alpha * B_p_k; 

	    		v_k = local_zeros(local_size(r));
	    		this->precond_->apply(r, v_k);

	    		g_v_prod_new = dot(r, v_k); 


	    		// if preconditioner yields nans or inf, or is precond. dir is indefinite - just return current step 
				if(!std::isfinite(g_v_prod_new)){
		    		return true; 
				}

	    		betta  = g_v_prod_new/ g_v_prod_old; 
	    		p_k = betta * p_k - v_k; 
 		
	    		// updating norms recursively  - see TR book
	    		sMp = (betta * sMp) + (alpha * p_norm); 
	    		p_norm = g_v_prod_new + (betta*betta * p_norm); 
	    		s_norm = s_norm_new; 


	    		g_norm = norm2(r);  

	    		if(this->verbose())
                    PrintInfo::print_iter_status(it, {g_norm, s_norm, p_norm, sMp});
	    		
	    		if(!std::isfinite(g_norm))
	    			return false; 

                converged = this->check_convergence(it, g_norm, 1, 1);
                it++;
	    	}

        	return true;
        }



	private:
	    void init(const SizeType &ls)
        {
            auto zero_expr = local_zeros(ls);

            //resets all buffers in case the size has changed
            if(!empty(v_k)) {
                v_k = zero_expr;
            }

            if(!empty(r)) {
                r = zero_expr;
            }

            if(!empty(p_k)) {
                p_k = zero_expr;
            }

            if(!empty(B_p_k)) {
                B_p_k = zero_expr;
            }
        }        


    private:
     	Vector v_k, r, p_k, B_p_k; 



    };


}

#endif //UTOPIA_TR_SUBPROBLEM_STEIHAUG_TOINT_HPP
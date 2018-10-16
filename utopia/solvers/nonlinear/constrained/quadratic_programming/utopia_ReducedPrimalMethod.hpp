#ifndef UTOPIA_REDUCED_PRIMAL_METHOD_HPP
#define UTOPIA_REDUCED_PRIMAL_METHOD_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Core.hpp"


#include <cmath>
#include <cassert>

namespace utopia 
{
	
	template<class Matrix, class Vector>
	class ReducedPrimalMethod : public IterativeSolver<Matrix, Vector> 
	{
	
	public:
		typedef utopia::BoxConstraints<Vector>  BoxConstraints;
		DEF_UTOPIA_SCALAR(Matrix)

		ReducedPrimalMethod()
		{
		}

		ReducedPrimalMethod(const ReducedPrimalMethod &) = default;

		inline ReducedPrimalMethod * clone() const override
		{
			return new ReducedPrimalMethod(*this);
		}		

		virtual bool set_box_constraints(const BoxConstraints & box)
		{
			constraints_ = box;
			constraints_.fill_empty_bounds();
			return true;
		}

		virtual void set_parameters(const Parameters params) override
		{
			IterativeSolver<Matrix, Vector>::set_parameters(params);
		}

		bool apply(const Vector &b, Vector &x) override
		{
			utopia_error("ReducedPrimalMethod:: apply not implemented yet.... \n"); 
			return false;
		}

		virtual void update(const std::shared_ptr<const Matrix> &op) override
		{
		    IterativeSolver<Matrix, Vector>::update(op);
		}


    Scalar compute_alpha_star(const Vector & x_cp, const Vector & lb, const Vector & ub, const Vector & d, const Vector & feasible_set) const
    {
        Vector alpha_stars = local_values(local_size(feasible_set).get(0), 1.0);  

        {
            Read<Vector>  rf(feasible_set); 
            Read<Vector>  rd(d); 
            Read<Vector>  ru(ub); 
            Read<Vector>  rlb(lb); 

            Write<Vector> wv(alpha_stars); 

            auto rr = range(alpha_stars); 

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                // TODO:: put approx eq
                if(feasible_set.get(i)==1)
                {
                    Scalar val = 1.0; 
                    Scalar di = d.get(i);  Scalar xi = x_cp.get(i); Scalar li = lb.get(i);  Scalar ui = ub.get(i); 

                    if(di!=0)
                        val = (di>0) ? (ui - xi)/di : (li - xi)/di; 

                    // checks for nans and infs
                    if(val < 1.0 && std::isfinite(val) && val > -1e9)
                        alpha_stars.set(i, val); 

                }
            }     
        }

        return min(alpha_stars); 
    }




    void prolongate_reduced_corr(const Vector & x_reduced,  const Vector & feasible_set,  Vector & x_prolongated) const
    {
        x_prolongated = local_values(local_size(feasible_set).get(0), 0.0); 


        {
            Write<Vector>  wx(x_prolongated); 

            Read<Vector>   rv(x_reduced); 
            Read<Vector>   rf(feasible_set); 

            auto rr = range(feasible_set); 
            auto rr_red = range(x_reduced); 

            SizeType counter = 0; 

            for(auto i = rr.begin(); i != rr.end(); ++i)
            {
                if(feasible_set.get(i)==1)
                {
                    Scalar val = x_reduced.get(rr_red.begin() + counter); 
                    x_prolongated.set(i, val); 
                    counter++; 
                }
            }
        }
    }

        void build_reduced_quantities(  const Vector &lb, const Vector & ub, const Vector & x_cp,
                                        const Vector & g, const Matrix & H, 
                                        Vector & feasible_set,  Matrix & H_reduced, Vector & g_reduced) const
        {
            feasible_set = local_values(local_size(lb).get(0), 0.); 

            // std::cout<<"build_reduced_quantities: 1 ---- \n"; 

            SizeType local_feasible_set = 0; 

            {
                Write<Vector>  rv(feasible_set); 

                Read<Vector>  rv1(ub); 
                Read<Vector>  rv2(lb); 
                Read<Vector>  rv3(x_cp); 

                auto rr = range(feasible_set); 

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {   
                    // TODO:: put approx eq
                    if(lb.get(i) != x_cp.get(i) &&  ub.get(i) != x_cp.get(i)){
                        feasible_set.set(i, 1); 
                        local_feasible_set++; 
                    }
                }       
            }

            // we can just get things directly 
            if(size(feasible_set).get(0)==sum(feasible_set))
            {
                g_reduced = g; 
                H_reduced = H; 

                return; 
            }


            g_reduced = local_zeros(local_feasible_set); 

            {
                SizeType local_counter = 0; 

                Read<Vector>  rf(feasible_set); 
                Read<Vector>  rg(g); 

                Write<Vector>  wg(g_reduced); 
                auto rr = range(feasible_set); 
                auto range_reduced = range(g_reduced); 

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {   
                    // TODO:: put approx eq
                    if(feasible_set.get(i)==1)
                    {
                        g_reduced.set(range_reduced.begin() + local_counter, g.get(i)); 
                        local_counter++;
                    }
                }       
            }

            // TODO:: use dense matrices for this kind of things 
            H_reduced  = local_values(local_size(H).get(0), local_feasible_set, 0.0); 

            if(local_size(feasible_set).get(0) != local_size(H).get(1))
                utopia_error("feasible_set, H: local sizes do not match .... \n"); 

            if(local_size(H_reduced).get(0) != local_size(H).get(0))
                utopia_error("H_reduced, H: local sizes do not match .... \n");             

            {
                Read<Vector>  rf(feasible_set); 
                Read<Matrix>  rM(H); 

                Write<Matrix>  wg(H_reduced); 

                auto row_original = row_range(H); 
                auto col_original = col_range(H); 

                // auto row_reduced = row_range(H_reduced); 
                auto col_reduced = col_range(H_reduced); 

                SizeType local_counter = 0; 

                // for (SizeType c = col_original.begin(); c != col_original.end(); ++c)
                for (SizeType r = row_original.begin(); r != row_original.end(); ++r)
                {   

                    for (SizeType c = col_original.begin(); c != col_original.end(); ++c)
                    {            
                        if(c==col_original.begin())
                            local_counter=0;

                        // std::cout<<"r: "<< r << "  c: "<< c << "  \n"; 

                        // TODO:: put approx eq
                        if(feasible_set.get(c)==1)
                        {
                            H_reduced.set(r, col_reduced.begin() + local_counter,  H.get(r, c)); 
                            local_counter++;
                        }
                    }    


                }  

            }
        }






	private:
		BoxConstraints constraints_;

	};
}

#endif //UTOPIA_REDUCED_PRIMAL_METHOD_HPP

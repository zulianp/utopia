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
		
		typedef UTOPIA_SCALAR(Vector)    Scalar;
    	typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////   CHECKED   /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        SizeType get_local_size_feasible_set(const Vector & feasible_set) const
        {
        	SizeType local_feasible_set = 0; 
            each_read(feasible_set, [&local_feasible_set](const SizeType i, const Scalar value) 
            {
				if(approxeq(value, 1.0)) 
					local_feasible_set++; 
			});

            return local_feasible_set; 
        }


        void build_feasible_set(const Vector &x, const Vector &ub, const Vector &lb, Vector & feasible_set) const 
        { 
        	if(empty(feasible_set) || local_size(feasible_set)!=local_size(x))
        		feasible_set = local_zeros(local_size(x)); 

			Read<Vector> r_ub(ub), r_lb(lb), r_x(x);
          	each_write(feasible_set, [ub, lb, x](const SizeType i) -> double { 
                      	Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  
                      	if(li < xi && xi < ui)
                        	return 1.0; 
                      	else
                        	return 0.0; }   );
        }



        void build_active_set(const Vector &x, const Vector &ub, const Vector &lb, Vector & active_set) const 
        { 
        	if(empty(active_set) || local_size(active_set)!=local_size(x))
        		active_set = local_zeros(local_size(x)); 

			Read<Vector> r_ub(ub), r_lb(lb), r_x(x);
          	each_write(active_set, [ub, lb, x](const SizeType i) -> double { 
                      	Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  
                      	if(li < xi && xi < ui)
                        	return 0.0; 
                      	else
                        	return 1.0; }   );
        }


        void build_reduced_vector(const Vector &x, const Vector & feasible_set, Vector &x_reduced) const
        {
        	if(local_size(x)==local_size(x_reduced) && size(x)==size(x_reduced))
        	{
        		x_reduced = x; 
        		return; 
        	}

            SizeType local_counter = 0; 
            {
	            Read<Vector>  rf(feasible_set); 
	            Read<Vector>  rg(x); 

	            Write<Vector>  wg(x_reduced); 
	            auto rr = range(feasible_set); 
	            auto range_reduced = range(x_reduced); 

	            for (SizeType i = rr.begin(); i != rr.end(); ++i)
	            {   
	                if(approxeq(feasible_set.get(i), 1.0))
	                {
	                    x_reduced.set(range_reduced.begin() + local_counter, x.get(i)); 
	                    local_counter++;
	                }
	            }   
	        }

        }



        void build_reduced_matrix(const Matrix &M, const Vector & feasible_set, Matrix &M_reduced) const 
        {
        	if(local_size(feasible_set).get(0) != local_size(M).get(0))
                utopia_error("feasible_set, H: local sizes do not match .... \n"); 

            if(local_size(M_reduced).get(1) != local_size(M).get(1))
                utopia_error("H_reduced, H: local sizes do not match .... \n");             


           	if(local_size(M)==local_size(M_reduced) && size(M)==size(M_reduced))
        	{
        		M_reduced = M; 
        		return; 
        	}


			{
                Read<Vector>  rf(feasible_set); 
                Read<Matrix>  rM(M); 

                Write<Matrix>  wg(M_reduced); 

                auto row_original = row_range(M); 
                auto col_original = col_range(M); 

                auto row_reduced = row_range(M_reduced); 

                SizeType local_row_counter = 0; 

                for (SizeType r = row_original.begin(); r != row_original.end(); ++r)
                {   
                	if(approxeq(feasible_set.get(r), 1.0))
                	{
                		for (SizeType c = col_original.begin(); c != col_original.end(); ++c)
	                    {            
	                        M_reduced.set(row_reduced.begin() + local_row_counter, c,  M.get(r, c)); 
	                    }    
                		local_row_counter++; 
                	}
                }  
            }


        }


	    void prolongate_reduced_corr(const Vector & x_reduced,  const Vector & feasible_set,  Vector & x_prolongated) const
	    {
	    	if(local_size(feasible_set) != local_size(x_prolongated) || size(feasible_set) != size(x_prolongated))
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


	private:
		BoxConstraints constraints_;

	};
}

#endif //UTOPIA_REDUCED_PRIMAL_METHOD_HPP

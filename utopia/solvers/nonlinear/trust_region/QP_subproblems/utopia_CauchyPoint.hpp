#ifndef UTOPIA_TR_SUBPROBLEM_CAUCHY_POINT_HPP
#define UTOPIA_TR_SUBPROBLEM_CAUCHY_POINT_HPP
#include "utopia_TRSubproblem.hpp"


namespace utopia 
{

	/**
	 * @brief   
           	This class evaluates the Cauchy Point. 
			\details 
				\f$ p_k^C = - \tau_k \frac{\Delta_k}{||g_k||} g_k \f$, where 
                \f$ \tau_k  =         \begin{cases}
                1 &      g_k^T B_k g_k \leq 0   \\
                \min(\frac{ ||g_k||^3 }{\Delta_k g_k^T B_k g_k}, 1 )& otherwise
            \end{cases}  \f$
	 */	
	template<class Matrix, class Vector>
    class CauchyPoint final: public TRSubproblem<Matrix, Vector>
    {
    	typedef UTOPIA_SCALAR(Vector) Scalar;

    public:

    	CauchyPoint(const Parameters params = Parameters()): TRSubproblem<Matrix, Vector>(params) {}

        bool apply(const Vector &b, Vector &x) override
        {
            return aux_solve(*this->get_operator(), -1.0 * b, x);
        }

        inline CauchyPoint * clone() const override
        {
            return new CauchyPoint();
        }        

    private:
        bool aux_solve(const Matrix &B, const Vector &g, Vector &p_k)
        {
			Scalar g_norm = norm2(g);
	    	Scalar g_B_g = dot(g, B * g); 
	    	Scalar tau = std::min(1.0, pow(g_norm,3.0) / (this->current_radius() * g_B_g));
	    	p_k = -tau * (this->current_radius() / g_norm) * (g) ;
	    	return true; 
        }

    };
}


#endif //UTOPIA_TR_SUBPROBLEM_CAUCHY_POINT_HPP


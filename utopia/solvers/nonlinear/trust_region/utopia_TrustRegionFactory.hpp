#ifndef UTOPIA_TR_STRATEGY_FACTORY_HPP
#define UTOPIA_TR_STRATEGY_FACTORY_HPP

#include "utopia_Core.hpp"
#include "utopia_Parameters.hpp"
#include "utopia_Function.hpp"
#include "utopia_FunctionNormalEq.hpp"
#include "utopia_TRNormalEquation.hpp"
#include "utopia_TrustRegion.hpp"

namespace utopia 
{

	typedef const char * TRStrategyTag;

	static TRStrategyTag AUTO_TR_TAG 			= "AUTO_TR";
	static TRStrategyTag CAUCHYPOINT_TAG 		= "CAUCHYPOINT";
	static TRStrategyTag DOGLEG_TAG 			= "DOGLEG";

	static TRStrategyTag STEIHAUG_TOINT_TAG 	= "STEIHAUG_TOINT";
	static TRStrategyTag TOINT_TAG 				= "TOINT";
	static TRStrategyTag NASH_TAG 		    	= "NASH";
	static TRStrategyTag LANCZOS_TAG 		    = "LANCZOS";
	static TRStrategyTag CGNE_TAG 		    	= "CGNE";


	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
	class TRStrategyFactory {
	public:
		TRStrategyFactory()
		{
			static_assert(Backend < HOMEMADE, "TRStrategyFactory not implemented for this backend");
		}
	};

	template<class Matrix, class Vector>
	typename TRStrategyFactory<Matrix, Vector>::StrategyPtr trust_region_strategy(const TRStrategyTag &tag = AUTO_TAG)
	{
	 	return TRStrategyFactory<Matrix, Vector>::new_trust_region_strategy(std::string(tag));
	}

	/**
	 * @brief      Trust region solve function
	 *
	 * @ingroup 	non-linear
	 * 
	 * 
	 * 				Note that TR solver can be used in 2 different settings: minimization of the objective function, solving nonlinear eq. \n
	 * 				Setting Function and relative derivatives properly determines which strategy to use. \n
	 * 				Using LeastSquaresFunction helps to solve arising LS in "smart" way, without explicitly forming hessian matrices with doubled conditioned number.  \n
	 *  			Currently, we have available: 
	 * 				-# <a href="http://www.mcs.anl.gov/~anitescu/CLASSES/2012/LECTURES/S310-2012-lect5.pdf">Cauchy point</a>
	 * 				-# <a href="http://www.numerical.rl.ac.uk/people/nimg/course/lectures/raphael/lectures/lec7slides.pdf">Dogleg</a>
	 * 				-# <a href="http://www.ccom.ucsd.edu/~peg/papers/trust.pdf">Steihaug-Toint</a> 	
	 * 				-# <a href="http://www.machinelearning.org/proceedings/icml2007/papers/114.pdf">Nash</a>
	 * 				-# <a href="http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.49.9338&rep=rep1&type=pdf">Generalized Lanczos</a>  
	 * 				-# <a href= "http://www.ii.uib.no/~trond/publications/papers/trust.pdf"> Toint</a> 
	 * 				-# <a href="\http://epubs.siam.org/doi/abs/10.1137/0613049?journalCode=sjmael">CGNE</a>, which is CG for normal eq.	
	 * @param      fun     The fun with nonlinear application context. 
	 * @param      x       Initial guess/solution
	 * @param[in]  params  The parameters
	 * 	 
	 */
	template<typename Matrix, typename Vector>
	const Parameters trust_region_solve(Function<Matrix, Vector> &fun, Vector &x, const Parameters params = Parameters())
	{
		// TODO:: check options proper options for normal eq. 
		if(LeastSquaresFunction<Matrix, Vector> * fun_ne_ptr = dynamic_cast<LeastSquaresFunction<Matrix, Vector>* >(&fun))
		{
			auto subproblem = trust_region_strategy<Matrix, Vector>(params.trust_region_alg()); 
			LeastSquaresTrustRegion<Matrix, Vector> tr_solver(subproblem);
			
			tr_solver.set_parameters(params);  
        	tr_solver.solve(*fun_ne_ptr, x);  
        	return tr_solver.parameters(); 
		}
		else
		{
			auto subproblem = trust_region_strategy<Matrix, Vector>(params.trust_region_alg()); 
			TrustRegion<Matrix, Vector> tr_solver(subproblem);
			
			tr_solver.set_parameters(params);  
        	tr_solver.solve(fun, x);  
        	return tr_solver.parameters(); 
		}

	}

}


#endif //UTOPIA_TR_STRATEGY_FACTORY_HPP

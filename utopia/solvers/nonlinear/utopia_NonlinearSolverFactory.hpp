#ifndef UTOPIA_NONLINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_NONLINEAR_SOLVER_FACTORY_HPP

#include "utopia_Core.hpp"
#include "utopia_FactoryMethod.hpp"

namespace utopia {
	typedef const char * NonlinearSolverTag;

	static NonlinearSolverTag NEWTON_TAG 			= "NEWTON";
	static NonlinearSolverTag LINE_SEARCH_TAG 		= "LINE_SEARCH";
	static NonlinearSolverTag TRUST_REGION_TAG 		= "TRUST_REGION";

	////////////////////////////////// depreciated /////////////////////////////////////
	/**
	 * @brief      Front-end to create nonlinear solver objects.
	 *
	 * @tparam     Matrix  
	 * @tparam     Vector  
	 */
	template<typename Matrix, typename Vector>
	class NonlinearSolverFactory {
	public: 

		typedef utopia::NewtonBase<Matrix, Vector> NonLinearSolverT;
		typedef std::shared_ptr<NonLinearSolverT> NonLinearSolverPtr; 
		typedef utopia::IFactoryMethod<NonLinearSolverT> FactoryMethodT;

		template<class Alg>
		using NLSolverFactoryMethod = FactoryMethod<NonLinearSolverT, Alg>;

		typedef std::shared_ptr< LinearSolver<Matrix, Vector> >  LinearSolverPtr;
		std::map<std::string, FactoryMethodT> nl_solvers_;

		inline static NonLinearSolverPtr solver(const NonlinearSolverTag &tag, LinearSolverPtr & linear_solver)
		{
			auto it = instance().nl_solvers_.find(tag);
			if(it == instance().nl_solvers_.end()) {
				std::cout<<"LinearSolver not available, solving with Newton.  \n";
				return std::make_shared<Newton<Matrix, Vector> >(linear_solver);  
			} else {
				auto ptr = it->second->make();
				ptr->set_linear_solver(linear_solver);
				return ptr;
			}
		}

	private:
		inline static const NonlinearSolverFactory &instance()
		{
			static NonlinearSolverFactory instance_;
			return instance_;
		}

		NonlinearSolverFactory()
		{
			init();
		}

		void init()
		{
			nl_solvers_[NEWTON_TAG] 	  = std::make_shared< NLSolverFactoryMethod<Newton<Matrix, Vector>> >();   
			nl_solvers_[TRUST_REGION_TAG] = std::make_shared< NLSolverFactoryMethod<TrustRegionBase<Matrix, Vector>> >();  
		}
	};

////////////////////////////////////////////////////////////////////////////////////////////////////////////



	/** \addtogroup non-linear
	 * @brief Solve functions for non-linear systems
	 * @ingroup solving
	 * 
	 */

	/**
	 * @brief      	Function to solve nonlinear system. 
	 *
	 *
	 * @ingroup 	non-linear	
	 * @param      fun     The function with nonlinear application context. 
	 * @param      x       The initial gues/solution. 
	 * @param[in]  params  The parameters.
	 *
	 */
	template<class Matrix, class Vector>
	const Parameters solve(Function<Matrix, Vector> &fun, Vector &x, const Parameters params = Parameters())
	{
		if(params.solver_type() == TRUST_REGION_TAG)
			return trust_region_solve(fun, x, params); 
		else if(params.solver_type() == LINE_SEARCH_TAG)
			return line_search_solve(fun, x, params); 
		else
			return newton_solve(fun, x, params); 
	}


	/**
	 * @brief      	Newton solver/Newton solve with dumping parameter. 
	 * @ingroup 	non-linear
	 * @param      fun     The fun
	 * @param      x       The initial guess
	 * @param[in]  params  The parameters
	 *
	 * @return     Parameters containing convergence history.
	 */
	template<class Matrix, class Vector>
	const Parameters newton_solve(Function<Matrix, Vector> &fun, Vector &x, const Parameters params = Parameters())
	{
		auto lin_solver = linear_solver<Matrix, Vector>(params.lin_solver_type());
		lin_solver->set_parameters(params); 

		Newton<Matrix, Vector> nlsolver(lin_solver);
		
		nlsolver.set_parameters(params);  
		nlsolver.solve(fun, x);  
		return nlsolver.parameters(); 

	}



}


#endif //UTOPIA_NONLINEAR_SOLVER_FACTORY_HPP

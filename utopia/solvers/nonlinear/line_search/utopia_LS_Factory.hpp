// /*
// * @Author: alenakopanicakova
// * @Date:   2016-06-10
// * @Last Modified by:   alenakopanicakova
// * @Last Modified time: 2016-10-11
// */

#ifndef UTOPIA_LS_STRATEGY_FACTORY_HPP
#define UTOPIA_LS_STRATEGY_FACTORY_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolverFactory.hpp"

#include <map>
#include <string>
#include <memory>
#include "utopia_Parameters.hpp"

namespace utopia {

	typedef const char * LSStrategyTag;

	static LSStrategyTag SIMPLE_BACKTRACKING_TAG 	= "SIMPLE_BACKTRACKING";
	static LSStrategyTag BACKTRACKING_TAG 			= "BACKTRACKING";


	/**
	 * @brief      Front-end to create ls strategy objects.
	 *
	 * @tparam     Matrix  
	 * @tparam     Vector  
	 */
	template<typename Matrix, typename Vector>
	class LSStrategyFactory
	{
		public: 
			typedef std::shared_ptr< LSStrategy<Matrix, Vector> > StrategyPtr;
			std::map<std::string, StrategyPtr> strategies_;

			inline static StrategyPtr new_line_search_strategy(const LSStrategyTag &tag)
			{
				auto it = instance().strategies_.find(tag);
				if(it == instance().strategies_.end()) 
				{
					std::cout<<"Strategy not available, solving with Backtracking.  \n";  
					return std::make_shared<utopia::Backtracking<Matrix, Vector> >(); 
				} 
				else 
				{
					return it->second;
				}
			}

		private:
			inline static const LSStrategyFactory &instance()
			{
				static LSStrategyFactory instance_;
				instance_.init();
				return instance_;

			}

			void init()
			{
					strategies_[SIMPLE_BACKTRACKING_TAG] 		= std::make_shared<utopia::SimpleBacktracking<Matrix, Vector> >(); 
					strategies_[BACKTRACKING_TAG] 				= std::make_shared<utopia::Backtracking<Matrix, Vector> >(); 
			}
	};


	template<class Matrix, class Vector>
	typename LSStrategyFactory<Matrix, Vector>::StrategyPtr line_search_strategy(const LSStrategyTag &tag = AUTO_TAG)
	{
	 	return LSStrategyFactory<Matrix, Vector>::new_line_search_strategy(tag);
	}


	/**
	 * @brief      Function to solve nonlinear system with line search Newton solver. 
	 *
	 * @param      fun     The fun with nonlinear context. 
	 * @param      x       The initial guess/ solution.
	 * @param[in]  params  The parameters
	 *
	 *
	 *			DEFAULT LINE SEARCH PARAMETERS: 
	 *     		line_search_alg() = "BACKTRACKING";  \n
     *     		ls_rho() = 0.5; \n
     *     		c1() = 1e-4; \n
     *     		c2() = 1e-8; \n
	 *     		n_line_search_iters() = 50; \n
	 *     		line_search_inner_verbose() = true; \n
	 *     		alpha_min() = 1e-7; \n
	 */
	template<typename Matrix, typename Vector>
	const Parameters line_search_solve(Function<Matrix, Vector> &fun, Vector &x, const Parameters params = Parameters())
	{
		auto lin_solver = LinearSolverFactory<Matrix, Vector>::new_linear_solver(params.lin_solver_type());
		Newton<Matrix, Vector> ls_solver(lin_solver);
		
		auto strategy = line_search_strategy<Matrix, Vector>(params.line_search_alg()); 
		strategy->set_parameters(params);
			
        ls_solver.set_line_search_strategy(strategy);
        ls_solver.set_parameters(params);
        ls_solver.solve(fun, x);  
        return ls_solver.parameters(); 
	}
}


#endif //UTOPIA_LS_STRATEGY_FACTORY_HPP

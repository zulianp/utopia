// /*
// * @Author: alenakopanicakova
// * @Date:   2016-06-10
// * @Last Modified by:   alenakopanicakova
// * @Last Modified time: 2016-10-11
// */

#ifndef UTOPIA_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_LINEAR_SOLVER_FACTORY_HPP

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_Parameters.hpp"

namespace utopia  {

	typedef const char * SolverTag;
	static SolverTag AUTO_TAG 		= "AUTO";

	static SolverTag UTOPIA_CG_TAG = "UTOPIA_CG";

	static SolverTag CG_TAG 		= "CG";
	static SolverTag BICGSTAB_TAG 	= "BICGSTAB";
	static SolverTag DIRECT_TAG 	= "DIRECT";
	static SolverTag KSP_TAG 		= "KSP";

#ifdef WITH_LAPACK
	static SolverTag LAPACK_TAG 	= "LAPACK";
#endif //WITH_LAPACK

#ifdef WITH_UMFPACK
	static SolverTag UMFPACK_TAG 	= "UMFPACK";
#endif //WITH_UMFPACK

	/**
	 * @brief      Front-end to create linear solver objects.
	 *
	 * @tparam     Matrix   
	 * @tparam     Vector   
	 * @tparam     Backend  
	 */
	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
	class LinearSolverFactory {
	public:
		LinearSolverFactory() {
			static_assert(Backend < HOMEMADE, "LinearSolverFactory not available for Backend");
		}
	};


	/**
	 * @brief      Returns linear solver based on tag.
	 *
	 * @param[in]  tag     The tag - takes info on choice of LS. 
	 *
	 * @tparam     Matrix 
	 * @tparam     Vector 
	 *
	 * @return    
	 */
	template<class Matrix, class Vector>
	typename LinearSolverFactory<Matrix, Vector>::LinearSolverPtr linear_solver(const SolverTag &tag = AUTO_TAG)
	{
	 	return LinearSolverFactory<Matrix, Vector>::new_linear_solver(tag);
	}

	/** \addtogroup Linear 
	 * @brief Solve functions for linear systems.
	 * @ingroup solving
	 *  @{
	 */


	/**
	 * @brief      Solves the linear system A * x = rhs. If no params are provided, solver is choosen and set-up by default. 
	 *
	 * @param[in]  A       The A.
	 * @param[in]  rhs     The right hand side.
	 * @param      x       The initial guess/ solution.
	 * @param[in]  params  The parameters.
	 *
	 * @tparam     Matrix  
	 * @tparam     Vector  
	 *
	 * @return    
	 */
	template<class Matrix, class Vector>
	bool solve(const Matrix A, const Vector rhs, Vector &x, const Parameters params = Parameters ())
	{
		auto solver = linear_solver<Matrix, Vector>(params.lin_solver_type()); 
		solver->set_parameters(params); 
		solver->solve(A, rhs, x); 

		return true; 
	}

	/** @}*/
}

#endif //UTOPIA_LINEAR_SOLVER_FACTORY_HPP

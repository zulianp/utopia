// /*
// * @Author: alenakopanicakova
// * @Date:   2016-06-10
// * @Last Modified by:   alenakopanicakova
// * @Last Modified time: 2016-10-11
// */

#ifndef UTOPIA_UTOPIA_LINEAR_SOLVER_FACTORY2_HPP
#define UTOPIA_UTOPIA_LINEAR_SOLVER_FACTORY2_HPP

#include "utopia_Core.hpp"

namespace utopia 
{

	typedef const char * SolverTag;
	static SolverTag AUTO_TAG 		= "AUTO";

	static SolverTag UTOPIA_CG_TAG = "UTOPIA_CG";

	static SolverTag CG_TAG 		= "CG";
	static SolverTag BICGSTAB_TAG 	= "BICGSTAB";
	static SolverTag DIRECT_TAG 	= "DIRECT";
	static SolverTag KSP_TAG 		= "KSP";

	static SolverTag LAPACK_TAG 	= "LAPACK";
	static SolverTag UMFPACK_TAG 	= "UMFPACK";
	

	/**
	 * @brief      Front-end to create linear solver objects.
	 *
	 * @tparam     Matrix   
	 * @tparam     Vector   
	 * @tparam     Backend  
	 */
	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
	class LinearSolverFactory;

// -------------------------------------------PETSC------------------------------------------
	template<typename Matrix, typename Vector>
	class LinearSolverFactory<Matrix, Vector, PETSC>
	{
		public: 
			typedef std::shared_ptr< LinearSolver<Matrix, Vector> >  LinearSolverPtr;
			std::map<std::string, LinearSolverPtr> solvers_;

			inline static LinearSolverPtr new_linear_solver(const SolverTag &tag)
			{
				auto it = instance().solvers_.find(tag);
				if(it == instance().solvers_.end()) 
				{
					// std::cout<<"LinearSolver not available, solving with utopia_CG.  \n";  //TODO:: fix tests and put back
					return std::make_shared<ConjugateGradient<Matrix, Vector> >();
				} 
				else 
				{
					return it->second;
				}
			}

		private:

			inline static const LinearSolverFactory &instance()
			{
				static LinearSolverFactory instance_;
				instance_.init();
				return instance_;

			}

			void init()
			{
				#ifdef WITH_PETSC
					solvers_[UTOPIA_CG_TAG] = std::make_shared< ConjugateGradient<Matrix, Vector, HOMEMADE> >();
					solvers_[AUTO_TAG] 		= std::make_shared< BiCGStab<Matrix, Vector> >();
					solvers_[CG_TAG] 		= std::make_shared< ConjugateGradient<Matrix, Vector> >();
					solvers_[BICGSTAB_TAG] 	= std::make_shared< BiCGStab<Matrix, Vector> >();
					solvers_[KSP_TAG] 		= std::make_shared< KSPSolver<Matrix, Vector> >();
					solvers_[DIRECT_TAG] 	= std::make_shared< Factorization<Matrix, Vector> >();
				#endif // WITH_PETSC
			}
	};


// -------------------------------------------BLAS------------------------------------------
	template<typename Matrix, typename Vector>
	class LinearSolverFactory<Matrix, Vector, BLAS>
	{

		public: 
			typedef std::shared_ptr< LinearSolver<Matrix, Vector> >  LinearSolverPtr;
			std::map<std::string, LinearSolverPtr> solvers_;

			inline static LinearSolverPtr new_linear_solver(const SolverTag &tag)
			{
				auto it = instance().solvers_.find(tag);
				if(it == instance().solvers_.end()) 
				{
					//std::cout<<"LinearSolver not avaialble, solve with CG.  \n";   // TODO fix tests and put back
					return std::make_shared<ConjugateGradient<Matrix, Vector> >();
				} 
				else 
				{
					return it->second;
				}
			}

		private:

			inline static const LinearSolverFactory &instance()
			{
				static LinearSolverFactory instance_;
				instance_.init();
				return instance_;

			}

			
			void init()
			{
#ifdef WITH_LAPACK
					//FIXME
					solvers_[LAPACK_TAG] = std::make_shared< LUDecomposition<Matrix, Vector> >();
					solvers_[AUTO_TAG]   = std::make_shared< LUDecomposition<Matrix, Vector> >();
#endif //WITH_LAPACK
	

//FIXME does not work for dense matrices				
				// #ifdef WITH_UMFPACK
				// 	solvers_[UMFPACK_TAG] = std::make_shared<UmfpackLU>();
				// #endif //WITH_UMFPACK
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


#endif //UTOPIA_UTOPIA_LINEAR_SOLVER_FACTORY2_HPP

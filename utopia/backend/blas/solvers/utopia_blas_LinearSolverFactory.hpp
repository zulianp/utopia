#ifndef UTOPIA_BLAS_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_BLAS_LINEAR_SOLVER_FACTORY_HPP 

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_ConjugateGradient.hpp"

#ifdef WITH_LAPACK
#include "utopia_Lapack.hpp"
#endif //WITH_LAPACK

#ifdef WITH_UMFPACK
#include "utopia_UmfpackLU.hpp"
#endif //WITH_UMFPACK

#include <map>
#include <string>
#include <memory>

namespace utopia {

	// -------------------------------------------BLAS------------------------------------------
		template<typename Matrix, typename Vector>
		class LinearSolverFactory<Matrix, Vector, BLAS> {
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
						return LinearSolverPtr(it->second->clone());
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
}

#endif //UTOPIA_BLAS_LINEAR_SOLVER_FACTORY_HPP

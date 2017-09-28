#ifndef UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP

#include "utopia_Traits.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_ConjugateGradient.hpp"

#include "utopia_PETScKSPSolver.hpp"
#include "utopia_PETScKSPSolvers.hpp"

#include "utopia_PETScFactorization.hpp"
#include "utopia_PETScFactorizations.hpp"

#include <map>
#include <string>
#include <memory>

namespace utopia {

	template<typename Matrix, typename Vector>
	class LinearSolverFactory<Matrix, Vector, PETSC> {
	public: 
		typedef std::shared_ptr< LinearSolver<Matrix, Vector> >  LinearSolverPtr;
		std::map<std::string, LinearSolverPtr> solvers_;

		inline static LinearSolverPtr new_linear_solver(const SolverTag &tag)
		{
			auto it = instance().solvers_.find(tag);
			if(it == instance().solvers_.end()) {
				// std::cout<<"LinearSolver not available, solving with utopia_CG.  \n";  //FIXME fix tests and put back
				return std::make_shared<ConjugateGradient<Matrix, Vector> >();
			} else  {
				return it->second;
			}
		}

		inline void add_solver(const SolverTag &tag, const LinearSolverPtr &solver_ptr)
		{
			solvers_[tag] = solver_ptr;
		}

		inline static LinearSolverFactory &instance()
		{
			static LinearSolverFactory instance_;
			instance_.init();
			return instance_;
		}

	private:
		void init()
		{
			add_solver(UTOPIA_CG_TAG,  std::make_shared< ConjugateGradient<Matrix, Vector, HOMEMADE> >());
			add_solver(AUTO_TAG,  	   std::make_shared< BiCGStab<Matrix, Vector> >());
			add_solver(CG_TAG,  	   std::make_shared< ConjugateGradient<Matrix, Vector> >());
			add_solver(BICGSTAB_TAG,   std::make_shared< BiCGStab<Matrix, Vector> >());
			add_solver(KSP_TAG,  	   std::make_shared< KSPSolver<Matrix, Vector> >());
			add_solver(DIRECT_TAG,     std::make_shared< Factorization<Matrix, Vector> >());
		}
	};
}

#endif //UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP

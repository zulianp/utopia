#ifndef UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP

#include "utopia_Traits.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_ConjugateGradient.hpp"

#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_KSPSolvers.hpp"

#include "utopia_petsc_Factorization.hpp"
#include "utopia_petsc_Factorizations.hpp"
#include "utopia_FactoryMethod.hpp"

#include <map>
#include <string>
#include <memory>

namespace utopia {

	template<typename Matrix, typename Vector>
	class LinearSolverFactory<Matrix, Vector, PETSC> {
	public: 
		typedef utopia::LinearSolver<Matrix, Vector> LinearSolverT;
		typedef std::shared_ptr<LinearSolverT>  LinearSolverPtr;
		typedef utopia::IFactoryMethod<LinearSolverT> FactoryMethodT;

		template<class Alg>
		using LSFactoryMethod = FactoryMethod<LinearSolverT, Alg>;
		std::map<std::string, std::shared_ptr<FactoryMethodT> > solvers_;

		inline static LinearSolverPtr new_linear_solver(const SolverTag &tag)
		{
			auto it = instance().solvers_.find(tag);
			if(it == instance().solvers_.end()) {
				return std::make_shared<ConjugateGradient<Matrix, Vector> >();
			} else  {
				return it->second->make();
			}
		}

		inline static LinearSolverFactory &instance()
		{
			static LinearSolverFactory instance_;
			return instance_;
		}

	private:
		LinearSolverFactory()
		{
			init();
		}
		
		void init()
		{
			solvers_[UTOPIA_CG_TAG] = std::make_shared< LSFactoryMethod<ConjugateGradient<Matrix, Vector, HOMEMADE>> >();
			solvers_[AUTO_TAG]   	= std::make_shared< LSFactoryMethod<BiCGStab<Matrix, Vector>> >();
			solvers_[CG_TAG]   	    = std::make_shared< LSFactoryMethod<ConjugateGradient<Matrix, Vector>> >();
			solvers_[BICGSTAB_TAG]  = std::make_shared< LSFactoryMethod<BiCGStab<Matrix, Vector>> >();
			solvers_[KSP_TAG] 		= std::make_shared< LSFactoryMethod<KSPSolver<Matrix, Vector>> >();
			solvers_[DIRECT_TAG]    = std::make_shared< LSFactoryMethod<Factorization<Matrix, Vector>> >();
		}
	};

}

#endif //UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP

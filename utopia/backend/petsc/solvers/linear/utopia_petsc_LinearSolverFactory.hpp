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
#include "utopia_petsc_GMRES.hpp"

#include "utopia_make_unique.hpp"

#include <map>
#include <string>
#include <memory>

namespace utopia {

    template<typename Matrix, typename Vector>
    class LinearSolverFactory<Matrix, Vector, PETSC> {
    public:
        typedef utopia::LinearSolver<Matrix, Vector>  LinearSolverT;
        typedef std::unique_ptr<LinearSolverT>        LinearSolverPtr;
        typedef utopia::IFactoryMethod<LinearSolverT> FactoryMethodT;

        template<class Alg>
        using LSFactoryMethod = FactoryMethod<LinearSolverT, Alg>;
        std::map<std::string, std::shared_ptr<FactoryMethodT> > solvers_;

        inline static LinearSolverPtr new_linear_solver(const std::string &tag)
        {
            auto it = instance().solvers_.find(tag);
            if(it == instance().solvers_.end()) {
                return utopia::make_unique<ConjugateGradient<Matrix, Vector> >();
            } else  {
                return it->second->make();
            }
        }

        inline static LinearSolverPtr default_linear_solver()
        {
            return utopia::make_unique<GMRES<Matrix, Vector>>("bjacobi");
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
            solvers_[Solver::utopia_cg()] = utopia::make_unique< LSFactoryMethod<ConjugateGradient<Matrix, Vector, HOMEMADE>> >();
            solvers_[Solver::automatic()]   	= utopia::make_unique< LSFactoryMethod<BiCGStab<Matrix, Vector>> >();
            solvers_[Solver::cg()]   	    = utopia::make_unique< LSFactoryMethod<ConjugateGradient<Matrix, Vector>> >();
            solvers_[Solver::bicgstab()]  = utopia::make_unique< LSFactoryMethod<BiCGStab<Matrix, Vector>> >();
            solvers_[Solver::ksp()] 		= utopia::make_unique< LSFactoryMethod<KSPSolver<Matrix, Vector>> >();
            solvers_[Solver::direct()]    = utopia::make_unique< LSFactoryMethod<Factorization<Matrix, Vector>> >();
            solvers_["gmres"]       = utopia::make_unique< LSFactoryMethod<GMRES<Matrix, Vector>> >();
        }
    };

}

#endif //UTOPIA_PETSC_LINEAR_SOLVER_FACTORY_HPP

#ifndef UTOPIA_BLAS_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_BLAS_LINEAR_SOLVER_FACTORY_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_FactoryMethod.hpp"

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
                typedef utopia::LinearSolver<Matrix, Vector> LinearSolverT;
                typedef std::shared_ptr<LinearSolverT>  LinearSolverPtr;
                typedef utopia::IFactoryMethod<LinearSolverT> FactoryMethodT;

                template<class Alg>
                using LSFactoryMethod = FactoryMethod<LinearSolverT, Alg>;
                std::map<std::string, std::shared_ptr<FactoryMethodT> > solvers_;


                inline static LinearSolverPtr new_linear_solver(const SolverType &tag)
                {
                    auto it = instance().solvers_.find(tag);
                    if(it == instance().solvers_.end()) {
                        return std::make_shared<ConjugateGradient<Matrix, Vector> >();
                    } else {
                        return LinearSolverPtr(it->second->make());
                    }
                }

            private:

                inline static const LinearSolverFactory &instance()
                {
                    static LinearSolverFactory instance_;
                    return instance_;
                }

                LinearSolverFactory()
                {
                    init();
                }

                void init()
                {
    #ifdef WITH_LAPACK
                        solvers_[Solver::direct()]    = std::make_shared< LSFactoryMethod< LUDecomposition<Matrix, Vector>> >();
                        solvers_[Solver::automatic()] = std::make_shared< LSFactoryMethod< LUDecomposition<Matrix, Vector>> >();
    #endif //WITH_LAPACK
                }
        };
}

#endif //UTOPIA_BLAS_LINEAR_SOLVER_FACTORY_HPP

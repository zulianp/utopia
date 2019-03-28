#ifndef UTOPIA_NONLINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_NONLINEAR_SOLVER_FACTORY_HPP

#include "utopia_Core.hpp"
#include "utopia_FactoryMethod.hpp"
#include "utopia_Newton.hpp"
#include "utopia_TrustRegion.hpp"
#include "utopia_SolverType.hpp"

namespace utopia {

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

        inline static NonLinearSolverPtr solver(const SolverType &tag, LinearSolverPtr & linear_solver)
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
            nl_solvers_[Solver::newton()] 	  = std::make_shared< NLSolverFactoryMethod<Newton<Matrix, Vector>> >();
            nl_solvers_[Solver::trust_region()] = std::make_shared< NLSolverFactoryMethod<TrustRegion<Matrix, Vector>> >();
        }
    };

////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * @brief      	Newton solver/Newton solve with dumping parameter.
     * @ingroup 	non-linear
     * @param      fun     The fun
     * @param      x       The initial guess
     * @param[in]  params  The parameters
     * @param[in]  params  Linear solver
     *
     * @return     Parameters containing convergence history.
     */
    template<class Matrix, class Vector>
    const SolutionStatus & newton_solve(Function<Matrix, Vector> &fun, Vector &x, Input & params, const std::shared_ptr <LinearSolver<Matrix, Vector> > &lin_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >())
    {
        Newton<Matrix, Vector> nlsolver(lin_solver);
        nlsolver.read(params);
        nlsolver.solve(fun, x);
        return nlsolver.solution_status();
    }


}


#endif //UTOPIA_NONLINEAR_SOLVER_FACTORY_HPP

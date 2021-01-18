#ifndef UTOPIA_LS_STRATEGY_FACTORY_HPP
#define UTOPIA_LS_STRATEGY_FACTORY_HPP

#include "utopia_Backtracking.hpp"
#include "utopia_Core.hpp"
#include "utopia_FactoryMethod.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_SimpleBacktracking.hpp"

#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace utopia {

    using SolverType = const char *;

    // static SolverType SIMPLE_BACKTRACKING_TAG = "SIMPLE_BACKTRACKING";
    // static SolverType BACKTRACKING_TAG = "BACKTRACKING";

    /**
     * @brief      Front-end to create ls strategy objects.
     *
     * @tparam     Vector
     */
    template <typename Vector>
    class LSStrategyFactory {
    public:
        using StrategyPtr = std::shared_ptr<LSStrategy<Vector>>;
        using FactoryMethodT = utopia::IFactoryMethod<LSStrategy<Vector>>;

        template <class Alg>
        using LSFactoryMethodT = FactoryMethod<LSStrategy<Vector>, Alg>;

        std::map<std::string, std::shared_ptr<FactoryMethodT>> strategies_;

        inline static StrategyPtr new_line_search_strategy(const SolverType &tag) {
            auto it = instance().strategies_.find(tag);
            if (it == instance().strategies_.end()) {
                utopia::out() << "Strategy not available, solving with Backtracking.  \n";
                return std::make_shared<utopia::Backtracking<Vector>>();
            } else {
                return it->second->make();
            }
        }

    private:
        inline static const LSStrategyFactory &instance() {
            static LSStrategyFactory instance_;
            instance_.init();
            return instance_;
        }

        void init() {
            strategies_[Solver::simple_backtracking()] =
                std::make_shared<LSFactoryMethodT<SimpleBacktracking<Vector>>>();
            strategies_[Solver::backtracking()] = std::make_shared<LSFactoryMethodT<Backtracking<Vector>>>();
        }
    };

    template <class Vector>
    typename LSStrategyFactory<Vector>::StrategyPtr line_search_strategy(const SolverType &tag = Solver::automatic()) {
        return LSStrategyFactory<Vector>::new_line_search_strategy(tag);
    }

    /**
     * @brief      Function to solve nonlinear system with line search Newton
     *solver.
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
    template <typename Matrix, typename Vector>
    const SolutionStatus &line_search_solve(Function<Matrix, Vector> &fun,
                                            Vector &x,
                                            const SolverType &tag,
                                            Input &params) {
        // auto lin_solver = LinearSolverFactory<Matrix,
        // Vector>::new_linear_solver(params.lin_solver_type());
        auto lin_solver = std::make_shared<ConjugateGradient<Matrix, Vector>>();
        Newton<Matrix, Vector> ls_solver(lin_solver);
        auto strategy = line_search_strategy<Vector>(tag);

        ls_solver.set_line_search_strategy(strategy);
        ls_solver.read(params);
        ls_solver.solve(fun, x);
        return ls_solver.solution_status();
    }
}  // namespace utopia

#endif  // UTOPIA_LS_STRATEGY_FACTORY_HPP

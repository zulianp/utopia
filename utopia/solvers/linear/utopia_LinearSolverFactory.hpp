#ifndef UTOPIA_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_LINEAR_SOLVER_FACTORY_HPP

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_SolverType.hpp"

namespace utopia {

/*
 * @brief      Front-end to create linear solver objects.
 *
 * @tparam     Matrix
 * @tparam     Vector
 * @tparam     Backend
 */
template <typename Matrix, typename Vector,
          int Backend = Traits<Matrix>::Backend>
class LinearSolverFactory {
 public:
  LinearSolverFactory() {
    static_assert(Backend < HOMEMADE,
                  "LinearSolverFactory not available for Backend");
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
template <class Matrix, class Vector>
typename LinearSolverFactory<Matrix, Vector>::LinearSolverPtr linear_solver(
    const SolverType &tag = Solver::automatic()) {
  return LinearSolverFactory<Matrix, Vector>::new_linear_solver(tag);
}

/** \addtogroup Linear
 * @brief Solve functions for linear systems.
 * @ingroup solving
 *  @{
 */

/**
 * @brief      Solves the linear system A * x = rhs. If no params are provided,
 * solver is choosen and set-up by default.
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
template <class Matrix, class Vector>
bool solve(const Matrix A, const Vector rhs, Vector &x,
           const SolverType &tag = Solver::automatic()) {
  auto solver = linear_solver<Matrix, Vector>(tag);
  solver->solve(A, rhs, x);

  return true;
}

/** @}*/
}  // namespace utopia

#endif  // UTOPIA_LINEAR_SOLVER_FACTORY_HPP

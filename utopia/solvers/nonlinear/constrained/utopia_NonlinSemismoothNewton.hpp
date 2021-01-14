/*! \file utopia_NonlinSemismoothNewton.hpp
    Nonlinear Semismooth Newton method
    Created by Alena Kopanicakova
*/

#ifndef UTOPIA_SOLVER_NONLINSEMISMOOTH_NEWTON_HPP
#define UTOPIA_SOLVER_NONLINSEMISMOOTH_NEWTON_HPP

#include <vector>
#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_NonLinearSolver.hpp"

namespace utopia {
/**
 * @brief      Nonlinear Semi-smooth solver.
 *
 * @tparam     Matrix
 * @tparam     Vector
 */
template <class Matrix, class Vector>
class NonlinSemismoothNewton : public NewtonBase<Matrix, Vector> {
  using Scalar = typename utopia::Traits<Vector>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;
  typedef typename NewtonBase<Matrix, Vector>::Solver Solver;
  using BoxConstraints = utopia::BoxConstraints<Vector>;

 public:
  NonlinSemismoothNewton(const std::shared_ptr<Solver> &linear_solver)
      : NewtonBase<Matrix, Vector>(linear_solver) {}

  bool solve(Function<Matrix, Vector> &fun, Vector &x_new) override {
    using namespace utopia;

    Scalar c = 1;
    SizeType iterations = 1;

    // const SizeType local_N = local_size(x_new).get(0);
    auto vec_layout = layout(x_new);
    auto mat_layout = square_matrix_layout(vec_layout);

    Vector lambda(vec_layout, 0.0);
    Vector Ginvg, d, g;
    Vector x_old = x_new;

    Matrix Ac, Ic, M;
    Matrix G;
    G.identity(mat_layout, 1.0);

    Matrix Hessian;
    Vector grad;

    Vector upbo;
    if (constraints_->has_upper_bound())
      upbo = *constraints_->upper_bound();
    else
      utopia::out() << "NonlinSemismoothNewton does not support other types at "
                       "the moment.... \n";

    this->linear_solve(G, upbo, Ginvg);

    bool converged = false;
    this->init_solver("NON-LINEAR - SEMISMOOTH NEWTON METHOD",
                      {" it. ", " err "});

    while (!converged) {
      // this is super expensive
      d = lambda + c * (G * x_new - upbo);

      //! active, inactive constraints
      if (is_sparse<Matrix>::value) {
        Ac.sparse(mat_layout, 1, 0);
        Ic.sparse(mat_layout, 1, 0);
      } else {
        Ac.dense(mat_layout, 0.0);
        Ic.dense(mat_layout, 0.0);
      }

      {
        Read<Vector> r(d);
        Write<Matrix> w_Ac(Ac);
        Write<Matrix> w_Ic(Ic);

        const Range rr = row_range(Ac);

        for (SizeType i = rr.begin(); i != rr.end(); i++) {
          if (d.get(i) > 0) {
            Ac.set(i, i, 1.0);
          } else {
            Ic.set(i, i, 1.0);
          }
        }
      }

      fun.gradient(x_new, grad);
      fun.hessian(x_new, Hessian);

      g = Hessian * x_new - grad;

      M = Ac + Ic * Hessian;
      this->linear_solve(M, (Ic * g + Ac * Ginvg), x_new);
      lambda = (g - Hessian * x_new);

      Scalar err = norm2(x_new - x_old);

      // print iteration status on every iteration
      if (this->verbose_) PrintInfo::print_iter_status(iterations, {err});

      // check convergence and print interation info
      converged = this->check_convergence(iterations, err, 1, 1);
      x_old = x_new;

      iterations++;
    }
    return true;
  }

  virtual bool set_box_constraints(const std::shared_ptr<BoxConstraints> &box) {
    constraints_ = box;
    return true;
  }

  virtual std::shared_ptr<BoxConstraints> get_box_constraints() const {
    return constraints_;
  }

 private:
  std::shared_ptr<BoxConstraints> constraints_;
};

}  // namespace utopia
#endif  // UTOPIA_SOLVER_NONLINSEMISMOOTH_NEWTON_HPP
#ifndef UTOPIA_JFNK_MULTIGRID_HPP
#define UTOPIA_JFNK_MULTIGRID_HPP
#include <utility>

#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelInterface.hpp"

namespace utopia {
/**
 * @brief      The class for Line-search multilevel optimization algorithm.
 *
 * @tparam     Matrix
 * @tparam     Vector
 */
template <class Matrix, class Vector>
class JFNK_Multigrid final
    : public NonlinearMultiLevelInterface<Matrix, Vector>,
      public OperatorBasedLinearSolver<Matrix, Vector> {
  using Scalar = typename utopia::Traits<Vector>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;
  using HessianApproximation = utopia::HessianApproximation<Vector>;
  using HessianApproxPtr = std::shared_ptr<HessianApproximation>;

  typedef utopia::ExtendedFunction<Matrix, Vector> Fun;
  using FunPtr = std::shared_ptr<Fun>;

  typedef utopia::MatrixFreeLinearSolver<Vector> LinSolver;
  using LinSolverPtr = std::shared_ptr<LinSolver>;

  typedef utopia::Transfer<Matrix, Vector> Transfer;

 public:
  JFNK_Multigrid(const SizeType &n_levels)
      : NonlinearMultiLevelInterface<Matrix, Vector>(n_levels),
        OperatorBasedLinearSolver<Matrix, Vector>() {}

  std::string name() override { return "JFNK_Multigrid"; }

  void init_memory() override { utopia::out() << "-------- to be done \n"; }

  JFNK_Multigrid *clone() const override {
    auto jfnk = new JFNK_Multigrid(this->n_levels());
    jfnk->set_smoother(this->mf_lin_solvers_[1]);
    jfnk->set_coarse_grid_solver(this->mf_lin_solvers_[0]);

    return jfnk;
  }

  void read(Input &in) override {
    NonlinearMultiLevelInterface<Matrix, Vector>::read(in);

    // if (_smoother) {
    //   in.get("smoother", *_smoother);
    // }
    // if (_coarse_solver) {
    //   in.get("coarse_solver", *_coarse_solver);
    // }
    // if (_ls_strategy) {
    //   in.get("ls_strategy", *_ls_strategy);
    // }
  }

  void print_usage(std::ostream &os) const override {
    NonlinearMultiLevelInterface<Matrix, Vector>::print_usage(os);

    // this->print_param_usage(os, "coarse_solver", "NewtonBase",
    //                         "Input parameters for coarse level QP
    // solvers.",
    //                         "-");
    // this->print_param_usage(os, "smoother", "NewtonBase",
    //                         "Input parameters for fine level QP solver.",
    //                         "-");
    // this->print_param_usage(os, "ls_strategy", "LSStrategy",
    //                         "Input parameters for line-search strategy.",
    //                         "-");
  }

  void update(const Operator<Vector> &op) override {
    const auto &layouts = this->local_level_layouts();
    memory_.init_memory(layouts);

    for (std::size_t l = 0; l != mf_lin_solvers_.size(); ++l) {
      assert(mf_lin_solvers_[l]);
      assert(hessian_approxs_[l]);
      mf_lin_solvers_[l]->init_memory(this->local_level_layouts_[l]);
      hessian_approxs_[l]->initialize(this->memory_.x[l], this->memory_.g[l]);
    }
  }

  bool solve(const Operator<Vector> &A, const Vector &rhs,
             Vector &sol) override {
    this->update(A);
    return this->apply(rhs, sol);
  }

  bool apply(const Vector &rhs, Vector &x) override {
    bool converged = false;
    SizeType it = 0, n_levels = this->n_levels();
    Scalar r_norm, r0_norm = 1, rel_norm = 1, energy;

    std::cout << "--------- yes, here, we  are  -------------- \n";
    exit(0);

    //     std::string header_message =
    //         this->name() + ": " + std::to_string(n_levels) + " levels";
    //     this->init_solver(header_message,
    //                       {" it. ", "|| grad ||", "r_norm", "Energy"});

    //     this->init_memory();

    //     Vector g(layout(x_h), 0.0);
    //     fine_fun.gradient(x_h, g);
    //     r0_norm = norm2(g);
    //     r_norm = r0_norm;

    //     fine_fun.value(x_h, energy);

    //     if (this->verbose())
    //       PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy});

    //     it++;

    //     while (!converged) {
    //       this->multiplicative_cycle(fine_fun, x_h, rhs, n_levels);

    // #ifdef CHECK_NUM_PRECISION_mode
    //       if (has_nan_or_inf(x_h) == 1) {
    //         x_h.zeros(layout(x_h));
    //         return true;
    //       }
    // #endif

    //       fine_fun.gradient(x_h, g);
    //       fine_fun.value(x_h, energy);

    //       r_norm = norm2(g);
    //       rel_norm = r_norm / r0_norm;

    //       // print iteration status on every iteration
    //       if (this->verbose())
    //         PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy});

    //       // check convergence and print interation info
    //       converged = this->check_convergence(it, r_norm, rel_norm, 1);
    //       it++;
    //     }

    //     this->print_statistics(it);

    // #ifdef CHECK_NUM_PRECISION_mode
    //     if (has_nan_or_inf(x_h) == 1) exit(0);
    // #endif

    return true;
  }

  bool set_coarse_grid_solver(const std::shared_ptr<LinSolver> &solver) {
    if (static_cast<SizeType>(mf_lin_solvers_.size()) != this->n_levels()) {
      mf_lin_solvers_.resize(this->n_levels());
    }

    mf_lin_solvers_[0] = solver;
    return true;
  }

  bool set_smoother(const std::shared_ptr<LinSolver> &solver) {
    if (static_cast<SizeType>(mf_lin_solvers_.size()) != this->n_levels()) {
      mf_lin_solvers_.resize(this->n_levels());
    }

    // starting from level 1 ....
    for (auto l = 1; l != this->n_levels(); ++l) {
      mf_lin_solvers_[l] = std::shared_ptr<LinSolver>(solver->clone());
    }

    return true;
  }

  bool set_functions(const std::vector<FunPtr> &level_functions) override {
    NonlinearMultiLevelInterface<Matrix, Vector>::set_functions(
        level_functions);

    hessian_approxs_.clear();

    for (auto l = 0; l < this->n_levels(); l++) {
      hessian_approxs_.push_back(
          std::make_shared<JFNK<Vector>>(this->function(l)));
    }

    return true;
  }

 private:
  bool multiplicative_cycle(Fun &fine_fun, Vector &u_l, const Vector &f,
                            const SizeType &l) {
    // Vector g_fine, g_restricted, g_coarse, u_2l, s_coarse, s_fine, u_init;
    // Scalar alpha;

    // this->make_iterate_feasible(fine_fun, u_l);

    // // PRE-SMOOTHING
    // smoothing(fine_fun, u_l, f, this->pre_smoothing_steps());
    // fine_fun.gradient(u_l, g_fine);

    // g_fine -= f;

    // this->transfer(l - 2).restrict(g_fine, g_restricted);
    // this->transfer(l - 2).project_down(u_l, u_2l);

    // this->make_iterate_feasible(this->function(l - 2), u_2l);
    // this->zero_correction_related_to_equality_constrain(this->function(l -
    // 2),
    //                                                     g_restricted);

    // this->function(l - 2).gradient(u_2l, g_coarse);

    // u_init = u_2l;
    // g_coarse -= g_restricted;  // tau correction - g_diff in rmtr

    // if (l == 2) {
    //   coarse_solve(this->function(0), u_2l, g_coarse);
    // } else {
    //   // recursive call into FAS - needs to be checked
    //   for (SizeType k = 0; k < this->mg_type(); k++) {
    //     SizeType l_new = l - 1;
    //     this->multiplicative_cycle(this->function(l - 2), u_2l, g_coarse,
    //                                l_new);
    //   }
    // }

    // s_coarse = u_2l - u_init;
    // this->transfer(l - 2).interpolate(s_coarse, s_fine);
    // this->zero_correction_related_to_equality_constrain(fine_fun, s_fine);

    // u_l += s_fine;

    // // POST-SMOOTHING
    // smoothing(fine_fun, u_l, f, this->post_smoothing_steps());

    return true;
  }

 private:
  std::vector<LinSolverPtr> mf_lin_solvers_;
  std::vector<HessianApproxPtr> hessian_approxs_;

  JFNKLevelMemory<Matrix, Vector> memory_;
};

}  // namespace utopia

#endif  // UTOPIA_JFNK_MULTIGRID_HPP

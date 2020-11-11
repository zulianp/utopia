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

  typedef utopia::OperatorBasedLinearSolver<Matrix, Vector> LinSolver;
  using LinSolverPtr = std::shared_ptr<LinSolver>;

  typedef utopia::Transfer<Matrix, Vector> Transfer;

 public:
  JFNK_Multigrid(const SizeType &n_levels)
      : NonlinearMultiLevelInterface<Matrix, Vector>(n_levels),
        OperatorBasedLinearSolver<Matrix, Vector>(),
        init_mem_(false) {}

  std::string name() override { return "JFNK_Multigrid"; }

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

  void update(const Vector &s, const Vector &y, const Vector &x,
              const Vector &g) {
    std::cout << "----- update ------- \n";

    hessian_approxs_[this->n_levels() - 1]->update(s, y, x, g);
    this->memory_.x[this->n_levels() - 1] = x;
    this->memory_.g[this->n_levels() - 1] = g;

    for (auto l = this->n_levels() - 1; l > 0; l--) {
      std::cout << "----- update ------- " << l << "          \n";
      this->transfer(l - 1).project_down(this->memory_.x[l],
                                         this->memory_.x[l - 1]);

      // this->make_iterate_feasible(this->function(l - 1),
      //                             this->memory_.x[l - 1]);

      // this->transfer(l - 1).restrict(this->memory_.g[l],
      //                                this->memory_.g[l - 1]);

      this->function(l - 1).gradient(this->memory_.x[l - 1],
                                     this->memory_.g[l - 1]);

      this->zero_correction_related_to_equality_constrain(this->function(l - 1),
                                                          memory_.g[l - 1]);

      hessian_approxs_[l - 1]->update(
          this->memory_.x[l - 1], this->memory_.g[l - 1],
          this->memory_.x[l - 1], this->memory_.g[l - 1]);
    }
  }

  void initialize(const Vector &x, const Vector &g) {
    if (init_mem_ == false) {
      init_memory();
    }

    this->memory_.x[this->n_levels() - 1] = x;
    this->memory_.g[this->n_levels() - 1] = g;

    hessian_approxs_[this->n_levels() - 1]->initialize(x, g);

    for (auto l = this->n_levels() - 1; l > 0; l--) {
      this->transfer(l - 1).project_down(this->memory_.x[l],
                                         this->memory_.x[l - 1]);

      // this->make_iterate_feasible(this->function(l - 1),
      //                             this->memory_.x[l - 1]);

      // this->transfer(l - 1).restrict(this->memory_.g[l],
      //                                this->memory_.g[l - 1]);

      this->function(l - 1).gradient(this->memory_.x[l - 1],
                                     this->memory_.g[l - 1]);

      this->zero_correction_related_to_equality_constrain(this->function(l - 1),
                                                          memory_.g[l - 1]);

      hessian_approxs_[l - 1]->initialize(this->memory_.x[l - 1],
                                          this->memory_.g[l - 1]);
    }
  }

  void update(const Operator<Vector> &op) override {}

  void init_memory() override {
    const auto &layouts = this->local_level_layouts();
    memory_.init_memory(layouts);

    for (std::size_t l = 0; l != mf_lin_solvers_.size(); ++l) {
      assert(mf_lin_solvers_[l]);
      assert(hessian_approxs_[l]);
      mf_lin_solvers_[l]->init_memory(this->local_level_layouts_[l]);
      // hessian_approxs_[l]->initialize(this->memory_.x[l],
      // this->memory_.g[l]);
    }

    init_mem_ = true;
  }

  bool solve(const Operator<Vector> &A, const Vector &rhs,
             Vector &sol) override {
    this->update(A);
    return this->apply(rhs, sol);
  }

  bool apply(const Vector &rhs, Vector &x) override {
    bool converged = false;
    SizeType it = 0, n_levels = this->n_levels();
    Scalar r_norm, r0_norm = 1, rel_norm = 1;

    std::string header_message =
        this->name() + ": " + std::to_string(n_levels) + " levels";
    this->init_solver(header_message, {" it. ", "|| grad ||", "r_norm"});

    //////////////////////////////////////////////////////////////////////
    ////////////////////////////  just some test ////////////////////////
    //////////////////////////////////////////////////////////////////////

    auto multiplication_action_finest =
        hessian_approxs_[n_levels - 1]->build_apply_H();

    multiplication_action_finest->apply(x, this->memory_.res[n_levels - 1]);
    this->memory_.res[n_levels - 1] -= rhs;

    r_norm = norm2(this->memory_.res[n_levels - 1]);
    r0_norm = r_norm;

    if (this->verbose()) PrintInfo::print_iter_status(it, {r_norm, rel_norm});

    it++;

    this->memory_.x[n_levels - 1] = x;
    this->memory_.rhs[n_levels - 1] = rhs;

    while (!converged) {
      this->multiplicative_cycle(n_levels - 1);

#ifdef CHECK_NUM_PRECISION_mode
      if (has_nan_or_inf(this->memory_.x[n_levels - 1]) == 1) {
        this->memory_.x[n_levels - 1].set(0.0);
        return true;
      }
#endif

      multiplication_action_finest->apply(this->memory_.x[n_levels - 1],
                                          this->memory_.res[n_levels - 1]);
      this->memory_.res[n_levels - 1] -= rhs;

      r_norm = norm2(this->memory_.res[n_levels - 1]);

      rel_norm = r_norm / r0_norm;

      // print iteration status on every iteration
      if (this->verbose()) PrintInfo::print_iter_status(it, {r_norm, rel_norm});

      // check convergence and print interation info
      converged = this->check_convergence(it, r_norm, rel_norm, 1);
      it++;
    }

    x = this->memory_.x[n_levels - 1];

    this->print_statistics(it);

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
  bool multiplicative_cycle(const SizeType &l) {
    // PRE-SMOOTHING
    this->level_solve(l, this->memory_.rhs[l], this->memory_.x[l],
                      this->pre_smoothing_steps());

    this->compute_residual(l, this->memory_.x[l]);

    this->transfer(l - 1).restrict(this->memory_.res[l],
                                   this->memory_.rhs[l - 1]);

    this->memory_.x[l - 1].set(0.0);
    if (l == 1) {
      this->zero_correction_related_to_equality_constrain(this->function(l - 1),
                                                          memory_.rhs[l - 1]);

      const SizeType coarse_grid_its = this->memory_.rhs[l - 1].size();

      this->level_solve(l - 1, this->memory_.rhs[l - 1], this->memory_.x[l - 1],
                        coarse_grid_its);
    } else {
      this->multiplicative_cycle(l - 1);
    }

    // interpolate
    this->transfer(l - 1).interpolate(this->memory_.x[l - 1],
                                      this->memory_.c[l]);

    // correct
    this->memory_.x[l] += this->memory_.c[l];

    // POST-SMOOTHING
    this->level_solve(l, this->memory_.rhs[l], this->memory_.x[l],
                      this->post_smoothing_steps());

    return true;
  }

  bool level_solve(const SizeType &level, const Vector &rhs, Vector &x,
                   const SizeType &sweeps = 1000) {
    auto multiplication_action = hessian_approxs_[level]->build_apply_H();
    mf_lin_solvers_[level]->verbose(false);
    mf_lin_solvers_[level]->max_it(sweeps);
    return mf_lin_solvers_[level]->solve(*multiplication_action, rhs, x);
  }

  void compute_residual(const SizeType &level, Vector &x) {
    auto multiplication_action = hessian_approxs_[level]->build_apply_H();

    multiplication_action->apply(x, this->memory_.res[level]);
    this->memory_.res[level] =
        this->memory_.rhs[level] - this->memory_.res[level];
  }

 private:
  std::vector<LinSolverPtr> mf_lin_solvers_;
  std::vector<HessianApproxPtr> hessian_approxs_;

  JFNKLevelMemory<Matrix, Vector> memory_;

  bool init_mem_;
};

}  // namespace utopia

#endif  // UTOPIA_JFNK_MULTIGRID_HPP

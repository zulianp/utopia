#ifndef UTOPIA_RMTR_BASE_HPP
#define UTOPIA_RMTR_BASE_HPP

#include "utopia_AuthoredWork.hpp"
#include "utopia_CiteUtopia.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"
#include "utopia_RMTRParams.hpp"

#include "utopia_LS_Strategy.hpp"
#include "utopia_Level.hpp"
#include "utopia_Linear.hpp"
#include "utopia_TRSubproblem.hpp"

#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_TRBase.hpp"

#include "utopia_FunEvalsIncludes.hpp"
#include "utopia_LevelMemory.hpp"

namespace utopia {

template <class Matrix, class Vector,
          MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
class RMTRBase : public NonlinearMultiLevelBase<Matrix, Vector>,
                 public RMTRParams<Vector>,
                 public AuthoredWork<Kopanicakova2020Recursive> {
 public:
  using Scalar = typename utopia::Traits<Vector>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;

  typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

  RMTRBase(const SizeType &n_levels)
      : NonlinearMultiLevelBase<Matrix, Vector>(n_levels),
        RMTRParams<Vector>(),
        AuthoredWork<Kopanicakova2020Recursive>(),
        ml_derivs_(n_levels),
        init_(false) {}

  ~RMTRBase() override = default;

  void read(Input &in) override {
    NonlinearMultiLevelBase<Matrix, Vector>::read(in);
    RMTRParams<Vector>::read(in);
  }

  void print_usage(std::ostream &os) const override {
    NonlinearMultiLevelBase<Matrix, Vector>::print_usage(os);
    RMTRParams<Vector>::print_usage(os);
  }

  virtual void reset() { init_ = false; }

 public:
  bool solve(Vector &x_h) override;

 protected:
  virtual bool multiplicative_cycle(const SizeType &level);
  virtual bool local_tr_solve(const SizeType &level,
                              const LocalSolveType &solve_type);

 protected:
  template <MultiLevelCoherence T = CONSISTENCY_LEVEL,
            enable_if_t<is_any<T, FIRST_ORDER, FIRST_ORDER_MGOPT, SECOND_ORDER,
                               GALERKIN>::value,
                        int> = 0>
  bool get_multilevel_hessian(const Fun &fun, const SizeType &level) {
    return ml_derivs_.compute_hessian(level, fun, memory_.x[level]);
  }

  template <MultiLevelCoherence T = CONSISTENCY_LEVEL,
            enable_if_t<is_any<T, FIRST_ORDER_DF, SECOND_ORDER_DF,
                               FIRST_ORDER_MULTIPLICATIVE_DF,
                               FIRST_ORDER_ADDITIVE_MULTIPLICATIVE_DF>::value,
                        int> = 0>
  bool get_multilevel_hessian(const Fun & /*fun*/, const SizeType & /*level*/) {
    return false;
  }

  bool get_multilevel_gradient(const Fun &fun, const SizeType &level,
                               const Vector &s_global) {
    return ml_derivs_.compute_gradient(level, fun, memory_.x[level], s_global);
  }

  bool get_multilevel_gradient(const Fun &fun, const SizeType &level,
                               const Vector &s_global,
                               const Scalar &energy_new_level_dep) {
    return ml_derivs_.compute_gradient(level, fun, memory_.x[level], s_global,
                                       energy_new_level_dep);
  }

  bool get_multilevel_gradient(const Fun &fun, const SizeType &level) {
    return ml_derivs_.compute_gradient(level, fun, memory_.x[level]);
  }

  Scalar get_multilevel_energy(const Fun &fun, const SizeType &level,
                               const Vector &s_global) {
    return ml_derivs_.compute_energy(level, fun, memory_.x[level], s_global);
  }

  Scalar get_multilevel_energy(const Fun &fun, const SizeType &level) {
    return ml_derivs_.compute_energy(level, fun, memory_.x[level]);
  }

  Scalar get_multilevel_gradient_energy(const Fun &fun, const SizeType &level,
                                        const Vector &s_global) {
    return ml_derivs_.compute_gradient_energy(level, fun, memory_.x[level],
                                              s_global);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <MultiLevelCoherence T = CONSISTENCY_LEVEL,
            enable_if_t<is_any<T, FIRST_ORDER, FIRST_ORDER_DF,
                               FIRST_ORDER_MGOPT>::value,
                        int> = 0>
  bool init_consistency_terms(const SizeType &level,
                              const Scalar &energy_fine_level_dep) {
    // UTOPIA_NO_ALLOC_BEGIN("RMTR::region111");
    // Restricted fine level gradient

    this->transfer(level - 1).restrict(this->ml_derivs_.g[level],
                                       this->ml_derivs_.g_diff[level - 1]);

    // Projecting current iterate to obtain initial iterate on coarser grid
    this->transfer(level - 1).project_down(this->memory_.x[level],
                                           this->memory_.x[level - 1]);
    // UTOPIA_NO_ALLOC_END();

    // UTOPIA_NO_ALLOC_BEGIN("RMTR::region112");
    if (!this->skip_BC_checks()) {
      this->make_iterate_feasible(this->function(level - 1),
                                  this->memory_.x[level - 1]);
    }
    // UTOPIA_NO_ALLOC_END();

    //----------------------------------------------------------------------------
    //    initializing coarse level (deltas, constraints, hessian approx, ...)
    //----------------------------------------------------------------------------
    UTOPIA_NO_ALLOC_BEGIN("RMTR::init_level");
    this->init_level(level - 1);
    UTOPIA_NO_ALLOC_END();

    //----------------------------------------------------------------------------
    //                  first order coarse level objective managment
    //----------------------------------------------------------------------------
    // UTOPIA_NO_ALLOC_BEGIN("RMTR::region114");
    this->function(level - 1).gradient(this->memory_.x[level - 1],
                                       this->ml_derivs_.g[level - 1]);
    // UTOPIA_NO_ALLOC_END();

    // UTOPIA_NO_ALLOC_BEGIN("RMTR::region1145");
    if (!this->skip_BC_checks()) {
      this->zero_correction_related_to_equality_constrain(
          this->function(level - 1), this->ml_derivs_.g_diff[level - 1]);
    }
    // UTOPIA_NO_ALLOC_END();

    // UTOPIA_NO_ALLOC_BEGIN("RMTR::region115");
    bool smoothness_flg = this->check_grad_smoothness()
                              ? this->recursion_termination_smoothness(
                                    this->ml_derivs_.g_diff[level - 1],
                                    this->ml_derivs_.g[level - 1], level - 1)
                              : true;

    this->memory_.help[level - 1] = this->ml_derivs_.g_diff[level - 1];
    this->ml_derivs_.g_diff[level - 1] -= this->ml_derivs_.g[level - 1];
    // UTOPIA_NO_ALLOC_END();

    // TODO:: verify:: NOT REALLY correct ....
    std::cout << "-------- fix energy evaluation ---------- \n";
    this->memory_.energy[level - 1] = energy_fine_level_dep;
    // storing for the first order consistency iteration
    this->ml_derivs_.g[level - 1] = this->memory_.help[level - 1];

    return smoothness_flg;
  }

  template <
      MultiLevelCoherence T = CONSISTENCY_LEVEL,
      enable_if_t<is_same<T, FIRST_ORDER_MULTIPLICATIVE_DF>::value, int> = 0>
  bool init_consistency_terms(const SizeType &level,
                              const Scalar &energy_fine_level_dep) {
    this->transfer(level - 1).restrict(this->ml_derivs_.g[level],
                                       this->ml_derivs_.g_diff[level - 1]);

    // Projecting current iterate to obtain initial iterate on coarser grid
    this->transfer(level - 1).project_down(this->memory_.x[level],
                                           this->memory_.x[level - 1]);

    if (!this->skip_BC_checks()) {
      this->make_iterate_feasible(this->function(level - 1),
                                  this->memory_.x[level - 1]);
    }

    //----------------------------------------------------------------------------
    //    initializing coarse level (deltas, constraints, hessian approx, ...)
    //----------------------------------------------------------------------------
    this->init_level(level - 1);

    //----------------------------------------------------------------------------
    //                  first order coarse level objective managment
    //----------------------------------------------------------------------------

    // todo:: investigate if correct
    // std::cout << "-------------- init const. terms 2 ------------- \n";
    this->compute_s_global(level, this->memory_.s_working[level]);

    // std::cout << "-------------- init const. terms 3 ------------- \n";
    // Scalar energy_fine = this->get_multilevel_energy(
    //     this->function(level), level, this->memory_.s_working[level]);
    Scalar energy_fine = energy_fine_level_dep;

    // std::cout << "-------------- init const. terms 4 ------------- \n";
    this->function(level - 1).gradient(this->memory_.x[level - 1],
                                       this->ml_derivs_.g[level - 1]);

    // std::cout << "-------------- init const. terms 5 ------------- \n";
    Scalar energy_coarse = 0.0;
    this->function(level - 1).value(this->memory_.x[level - 1], energy_coarse);

    // std::cout << "-------------- init const. terms 6 ------------- \n";

    // if (!this->skip_BC_checks()) {
    //   this->zero_correction_related_to_equality_constrain(
    //       this->function(level - 1), this->ml_derivs_.g_diff[level - 1]);
    // }

    // TODO:: pre-allocate
    Vector help1 = (1. / energy_coarse * this->ml_derivs_.g_diff[level - 1]);
    Scalar f_cc = (energy_fine / (energy_coarse * energy_coarse));
    Vector help2 = (f_cc * this->ml_derivs_.g[level - 1]);

    // self.g_diff[level] = (1./loss_coarse * grad_fine) -
    // (loss_fine/(loss_coarse*loss_coarse) * grad_coarse)

    // storing for the first order consistency iteration
    this->ml_derivs_.g[level - 1] = this->ml_derivs_.g_diff[level - 1];

    this->ml_derivs_.g_diff[level - 1] = help1 - help2;

    this->ml_derivs_.e_diff[level - 1] = energy_fine / energy_coarse;

    // as level is also zeroth order consistent
    this->memory_.energy[level - 1] = energy_fine;

    return true;
  }

  template <
      MultiLevelCoherence T = CONSISTENCY_LEVEL,
      enable_if_t<is_same<T, FIRST_ORDER_ADDITIVE_MULTIPLICATIVE_DF>::value,
                  int> = 0>
  bool init_consistency_terms(const SizeType &level,
                              const Scalar &energy_fine_level_dep) {
    this->transfer(level - 1).restrict(this->ml_derivs_.g[level],
                                       this->ml_derivs_.g_diff_m[level - 1]);

    this->ml_derivs_.g_diff_a[level - 1] = this->ml_derivs_.g_diff_m[level - 1];

    // Projecting current iterate to obtain initial iterate on coarser grid
    this->transfer(level - 1).project_down(this->memory_.x[level],
                                           this->memory_.x[level - 1]);

    if (!this->skip_BC_checks()) {
      this->make_iterate_feasible(this->function(level - 1),
                                  this->memory_.x[level - 1]);
    }

    if (!this->skip_BC_checks()) {
      this->zero_correction_related_to_equality_constrain(
          this->function(level - 1), this->ml_derivs_.g_diff_a[level - 1]);
      this->zero_correction_related_to_equality_constrain(
          this->function(level - 1), this->ml_derivs_.g_diff_m[level - 1]);
    }

    //----------------------------------------------------------------------------
    //    initializing coarse level (deltas, constraints, hessian approx, ...)
    //----------------------------------------------------------------------------
    this->init_level(level - 1);

    //----------------------------------------------------------------------------
    //                  first order coarse level objective managment
    //----------------------------------------------------------------------------

    // todo:: investigate if correct
    this->compute_s_global(level, this->memory_.s_working[level]);
    // Scalar energy_fine = this->get_multilevel_energy(
    //     this->function(level), level, this->memory_.s_working[level]);
    Scalar energy_fine = energy_fine_level_dep;

    this->function(level - 1).gradient(this->memory_.x[level - 1],
                                       this->ml_derivs_.g[level - 1]);
    Scalar energy_coarse = 0.0;
    this->function(level - 1).value(this->memory_.x[level - 1], energy_coarse);

    // TODO:: preallocate
    Vector help1 = (1. / energy_coarse * this->ml_derivs_.g_diff_m[level - 1]);
    Scalar f_cc = (energy_fine / (energy_coarse * energy_coarse));
    Vector help2 = (f_cc * this->ml_derivs_.g[level - 1]);

    // storing for the first order consistency iteration
    this->ml_derivs_.g[level - 1] = this->ml_derivs_.g_diff_m[level - 1];
    // true, as also 0th order consistent
    this->memory_.energy[level - 1] = energy_fine;

    this->ml_derivs_.g_diff_m[level - 1] = help1 - help2;
    this->ml_derivs_.e_diff_m[level - 1] = energy_fine / energy_coarse;

    this->ml_derivs_.g_diff_a[level - 1] -= this->ml_derivs_.g[level - 1];
    this->ml_derivs_.e_diff_m[level - 1] = energy_fine - energy_coarse;

    return true;
  }

  virtual void init_hess_app_terms(const SizeType & /* level */) {}

  template <MultiLevelCoherence T = CONSISTENCY_LEVEL,
            enable_if_t<is_same<T, SECOND_ORDER_DF>::value, int> = 0>
  bool init_consistency_terms(const SizeType &level) {
    // UTOPIA_NO_ALLOC_BEGIN("RMTR::region111");
    // Restricted fine level gradient
    this->transfer(level - 1).restrict(this->ml_derivs_.g[level],
                                       this->ml_derivs_.g_diff[level - 1]);

    // Projecting current iterate to obtain initial iterate on coarser grid
    this->transfer(level - 1).project_down(this->memory_.x[level],
                                           this->memory_.x[level - 1]);
    // UTOPIA_NO_ALLOC_END();

    // UTOPIA_NO_ALLOC_BEGIN("RMTR::region112");
    if (!this->skip_BC_checks()) {
      this->make_iterate_feasible(this->function(level - 1),
                                  this->memory_.x[level - 1]);
    }
    // UTOPIA_NO_ALLOC_END();

    //----------------------------------------------------------------------------
    //    initializing coarse level (deltas, constraints, hessian approx, ...)
    //----------------------------------------------------------------------------
    UTOPIA_NO_ALLOC_BEGIN("RMTR::init_level");
    this->init_level(level - 1);
    UTOPIA_NO_ALLOC_END();

    //----------------------------------------------------------------------------
    //                  first order coarse level objective managment
    //----------------------------------------------------------------------------
    // UTOPIA_NO_ALLOC_BEGIN("RMTR::region114");
    this->function(level - 1).gradient(this->memory_.x[level - 1],
                                       this->ml_derivs_.g[level - 1]);
    // UTOPIA_NO_ALLOC_END();

    // UTOPIA_NO_ALLOC_BEGIN("RMTR::region1145");
    if (!this->skip_BC_checks()) {
      this->zero_correction_related_to_equality_constrain(
          this->function(level - 1), this->ml_derivs_.g_diff[level - 1]);
    }
    // UTOPIA_NO_ALLOC_END();

    // UTOPIA_NO_ALLOC_BEGIN("RMTR::region115");
    bool smoothness_flg = this->check_grad_smoothness()
                              ? this->recursion_termination_smoothness(
                                    this->ml_derivs_.g_diff[level - 1],
                                    this->ml_derivs_.g[level - 1], level - 1)
                              : true;
    this->ml_derivs_.g_diff[level - 1] -= this->ml_derivs_.g[level - 1];
    // UTOPIA_NO_ALLOC_END();

    this->init_hess_app_terms(level);

    return smoothness_flg;
  }

  template <MultiLevelCoherence T = CONSISTENCY_LEVEL,
            enable_if_t<is_same<T, SECOND_ORDER>::value, int> = 0>
  bool init_consistency_terms(const SizeType &level) {
    // Restricted fine level gradient
    this->transfer(level - 1).restrict(this->ml_derivs_.g[level],
                                       this->ml_derivs_.g_diff[level - 1]);

    // Projecting current iterate to obtain initial iterate on coarser grid
    this->transfer(level - 1).project_down(this->memory_.x[level],
                                           this->memory_.x[level - 1]);

    if (!this->skip_BC_checks()) {
      this->make_iterate_feasible(this->function(level - 1),
                                  this->memory_.x[level - 1]);
    }

    //----------------------------------------------------------------------------
    //    initializing coarse level (deltas, constraints, hessian approx, ...)
    //----------------------------------------------------------------------------
    this->init_level(level - 1);

    //----------------------------------------------------------------------------
    //                  first order coarse level objective managment
    //----------------------------------------------------------------------------
    this->function(level - 1).gradient(this->memory_.x[level - 1],
                                       this->ml_derivs_.g[level - 1]);

    if (!this->skip_BC_checks()) {
      this->zero_correction_related_to_equality_constrain(
          this->function(level - 1), this->ml_derivs_.g_diff[level - 1]);
    }

    bool smoothness_flg = this->check_grad_smoothness()
                              ? this->recursion_termination_smoothness(
                                    this->ml_derivs_.g_diff[level - 1],
                                    this->ml_derivs_.g[level - 1], level - 1)
                              : true;
    this->ml_derivs_.g_diff[level - 1] -= this->ml_derivs_.g[level - 1];

    //----------------------------------------------------------------------------
    //                   second order coarse level objective managment
    //----------------------------------------------------------------------------
    // testing...
    // if(!this->deltaH_lagging() || level==this->n_levels()-1){
    if (!this->deltaH_lagging()) {
      this->get_multilevel_hessian(this->function(level), level);
    }
    this->transfer(level - 1).restrict(this->ml_derivs_.H[level],
                                       this->ml_derivs_.H_diff[level - 1]);

    if (!this->skip_BC_checks()) {
      this->zero_correction_related_to_equality_constrain_mat(
          this->function(level - 1), this->ml_derivs_.H_diff[level - 1]);
    }

    this->function(level - 1).hessian(this->memory_.x[level - 1],
                                      this->ml_derivs_.H[level - 1]);
    this->ml_derivs_.H_diff[level - 1] -= this->ml_derivs_.H[level - 1];

    std::cout << "-------------- add grad, energy storage -------- \n";

    return smoothness_flg;
  }

  template <MultiLevelCoherence T = CONSISTENCY_LEVEL,
            enable_if_t<is_same<T, GALERKIN>::value, int> = 0>
  bool init_consistency_terms(const SizeType &level) {
    UTOPIA_NO_ALLOC_BEGIN("RMTR::init_consistency_terms0");

    this->truncate(level);

    // Restricted fine level gradient
    this->transfer(level - 1).restrict(this->ml_derivs_.g[level],
                                       this->ml_derivs_.g_diff[level - 1]);

    if (!this->skip_BC_checks()) {
      this->zero_correction_related_to_equality_constrain(
          this->function(level - 1), this->ml_derivs_.g_diff[level - 1]);
    }

    // Projecting current iterate to obtain initial iterate on coarser grid
    this->transfer(level - 1).project_down(this->memory_.x[level],
                                           this->memory_.x[level - 1]);

    //----------------------------------------------------------------------------
    //    initializing coarse level (deltas, constraints, hessian approx, ...)
    //----------------------------------------------------------------------------
    this->init_level(level - 1);
    UTOPIA_NO_ALLOC_END();

    UTOPIA_NO_ALLOC_BEGIN("RMTR::hessian_comp2");
    this->get_multilevel_hessian(this->function(level), level);
    UTOPIA_NO_ALLOC_END();

    this->transfer(level - 1).restrict(this->ml_derivs_.H[level],
                                       this->ml_derivs_.H_diff[level - 1]);

    UTOPIA_NO_ALLOC_BEGIN("RMTR::init_consistency_terms1");
    if (!this->skip_BC_checks()) {
      this->zero_correction_related_to_equality_constrain_mat(
          this->function(level - 1), this->ml_derivs_.H_diff[level - 1]);
    }
    UTOPIA_NO_ALLOC_END();

    std::cout << "-------------- add grad, energy storage -------- \n";

    return true;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <MultiLevelCoherence T = CONSISTENCY_LEVEL,
            enable_if_t<is_any<T, FIRST_ORDER, FIRST_ORDER_DF,
                               FIRST_ORDER_MGOPT, SECOND_ORDER_DF>::value,
                        int> = 0>
  bool init_deriv_loc_solve(const Fun & /*fun*/, const SizeType &level,
                            const LocalSolveType &solve_type) {
    if (!(solve_type == PRE_SMOOTHING && level == this->n_levels() - 1)) {
      if (solve_type == PRE_SMOOTHING || solve_type == COARSE_SOLVE) {
        this->ml_derivs_.g[level] += this->ml_derivs_.g_diff[level];
        this->memory_.gnorm[level] = this->criticality_measure(level);
      }
    }

    return true;
  }

  template <
      MultiLevelCoherence T = CONSISTENCY_LEVEL,
      enable_if_t<is_same<T, FIRST_ORDER_MULTIPLICATIVE_DF>::value, int> = 0>
  bool init_deriv_loc_solve(const Fun &fun, const SizeType &level,
                            const LocalSolveType &solve_type) {
    if (!(solve_type == PRE_SMOOTHING && level == this->n_levels() - 1)) {
      if (solve_type == PRE_SMOOTHING || solve_type == COARSE_SOLVE) {
        // this->get_multilevel_gradient(fun, level);
        // std::cout << "-------------- beg  here ----------- \n";
        this->memory_.gnorm[level] = this->criticality_measure(level);
      }
    }

    return true;
  }

  template <
      MultiLevelCoherence T = CONSISTENCY_LEVEL,
      enable_if_t<is_same<T, FIRST_ORDER_ADDITIVE_MULTIPLICATIVE_DF>::value,
                  int> = 0>
  bool init_deriv_loc_solve(const Fun &fun, const SizeType &level,
                            const LocalSolveType &solve_type) {
    if (!(solve_type == PRE_SMOOTHING && level == this->n_levels() - 1)) {
      if (solve_type == PRE_SMOOTHING || solve_type == COARSE_SOLVE) {
        std::cout << "-------------- fix -------- \n";
        this->get_multilevel_gradient(fun, level);
        this->memory_.gnorm[level] = this->criticality_measure(level);
      }
    }

    return true;
  }

  template <MultiLevelCoherence T = CONSISTENCY_LEVEL,
            enable_if_t<is_same<T, SECOND_ORDER>::value, int> = 0>
  bool init_deriv_loc_solve(const Fun & /*fun*/, const SizeType &level,
                            const LocalSolveType &solve_type) {
    bool make_hess_updates = true;

    if (!(solve_type == PRE_SMOOTHING && level == this->n_levels() - 1)) {
      if (solve_type == PRE_SMOOTHING || solve_type == COARSE_SOLVE) {
        std::cout << "-------------- fix -------- \n";
        this->ml_derivs_.g[level] += this->ml_derivs_.g_diff[level];
        this->memory_.gnorm[level] = this->criticality_measure(level);

        this->ml_derivs_.H[level] += this->ml_derivs_.H_diff[level];
        make_hess_updates = false;
      }
    }

    return make_hess_updates;
  }

  template <MultiLevelCoherence T = CONSISTENCY_LEVEL,
            enable_if_t<is_same<T, GALERKIN>::value, int> = 0>
  bool init_deriv_loc_solve(const Fun & /*fun*/, const SizeType &level,
                            const LocalSolveType &solve_type) {
    std::cout << "-------------- fix -------- \n";
    if (!(solve_type == PRE_SMOOTHING && level == this->n_levels() - 1)) {
      if ((solve_type == PRE_SMOOTHING && level < this->n_levels() - 1) ||
          (solve_type == COARSE_SOLVE)) {
        this->ml_derivs_.g[level] = this->ml_derivs_.g_diff[level];
        this->memory_.gnorm[level] = this->criticality_measure(level);
      }
    }

    return true;
  }

  ///////////////////////////////////////////////////////////// HELPERS
  //////////////////////////////////////////////////////////////////////////////////////////////////

  void print_level_info(const SizeType &level) {
    if (this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE &&
        mpi_world_rank() == 0) {
      if (level == 0) {
        std::cout << this->yellow_;
        std::string solver_type = "COARSE SOLVE:: " + std::to_string(level);
        this->print_init_message(
            solver_type, {" it. ", "|| g ||", "   E + <g_diff, x>", "ared   ",
                          "  pred  ", "  rho  ", "  delta "});
      } else {
        std::cout << this->green_;
        std::string solver_type = "SMOOTHER:  " + std::to_string(level);
        this->print_init_message(
            solver_type, {" it. ", "|| g ||", "   E + <g_diff, x>", "ared   ",
                          "  pred  ", "  rho  ", "  delta "});
      }
    }
  }

  virtual void truncate(const SizeType & /*level*/) {}

  virtual void make_ml_iterate_feasible(const SizeType & /*level*/) {}

  /**
   * @brief      "Heuristics", which decides if hessian needs to be updated or
   * now
   *
   * @param[in]  g_new   new gradient
   * @param[in]  g_old   Old gradient
   * @param[in]  s       The correction
   * @param[in]  H       The hessian
   * @param[in]  rho     The rho
   * @param[in]  g_norm  Norm of gradient
   *
   */
  virtual bool update_hessian(const Vector &g_new, const Vector &g_old,
                              const Vector &s, const Matrix &H,
                              const Scalar &rho, const Scalar &g_norm) {
    // iteration is not sucessful enough
    if (rho > 0 && rho < this->hessian_update_eta()) return true;

    // get rid of allocations
    Vector help = g_new - g_old - H * s;

    // Hessian approx is relativelly poor
    return (norm2(help) > this->hessian_update_delta() * g_norm) ? true : false;
  }

  virtual bool update_level(const SizeType & /*level*/,
                            const Scalar & /*energy_new_level_dep*/) {
    return false;
  }

  virtual Scalar update_local_grad(const bool &make_grad_updates,
                                   const SizeType &level, Scalar &rho,
                                   Scalar &energy_new) {
    if (make_grad_updates) {
      // std::cout<<"grad updated... \n"; s
      // Vector g_old = memory_.g[level];
      UTOPIA_NO_ALLOC_BEGIN("RMTR::grad_computation");
      this->get_multilevel_gradient(this->function(level), level,
                                    this->memory_.s_working[level], energy_new);
      this->memory_.gnorm[level] = this->criticality_measure(level);
      UTOPIA_NO_ALLOC_END();

      // this is just a safety check
      if (!std::isfinite(this->memory_.gnorm[level])) {
        UTOPIA_NO_ALLOC_BEGIN("RMTR::region8");
        rho = 0;
        this->memory_.x[level] -=
            this->memory_.s[level];  // return iterate into its initial state
        this->compute_s_global(level, this->memory_.s_working[level]);
        this->get_multilevel_gradient(this->function(level), level,
                                      this->memory_.s_working[level]);
        this->memory_.gnorm[level] = this->criticality_measure(level);
        UTOPIA_NO_ALLOC_END();
      }

      // make_hess_updates =   this->update_hessian(memory_.g[level], g_old, s,
      // H, rho, g_norm);
    }

    return rho;
  }

  void compute_s_global(const SizeType &level, Vector &s_global) {
    if (empty(this->memory_.x_0[level])) {
      utopia_error("this should not happen, remove when done testing... \n");
    } else if (level < this->n_levels() - 1) {
      s_global = this->memory_.x[level] - this->memory_.x_0[level];
    }
    // we do not need to compute s_working on the fines level
  }

  ///////////////////////////////////////////// INITIALIZATIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // to overwrite
  void init_memory() override {
    const auto &layouts = this->local_level_layouts();
    const std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector>>> &funs =
        this->level_functions();

    ml_derivs_.init_memory(layouts, funs);
    memory_.init_memory(layouts);

    // init deltas to some default value...
    for (Scalar l = 0; l < this->n_levels(); l++) {
      this->memory_.delta[l] = this->delta0();
    }
  }

  void handle_equality_constraints() {
    const SizeType L = this->n_levels();

    for (SizeType i = 0; i < L - 1; ++i) {
      const auto &flags = this->function(i + 1).get_eq_constrains_flg();
      assert(!empty(flags));
      this->transfer(i).handle_equality_constraints(flags);
    }
  }

  virtual void init_level(const SizeType &level) {
    this->memory_.delta[level] = this->memory_.delta[level + 1];
  }

  virtual bool check_initialization() {
    if (static_cast<SizeType>(this->level_functions_.size()) !=
        this->n_levels()) {
      utopia_error(
          "utopia::RMTR:: number of level Functions and levels not equal. \n");
      return false;
    }
    if (static_cast<SizeType>(this->transfers_.size()) + 1 !=
        this->n_levels()) {
      utopia_error(
          "utopia::RMTR:: number of transfers and levels not equal. \n");
      return false;
    }

    return true;
  }

  virtual bool check_feasibility(const SizeType & /*level */) { return false; }

  ///////////////////////////////////////////// CONVERGENCE CHECKS
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool check_global_convergence(const SizeType &it, const Scalar &r_norm,
                                const Scalar &rel_norm, const Scalar &delta) {
    bool converged = NonlinearMultiLevelBase<Matrix, Vector>::check_convergence(
        it, r_norm, rel_norm, 1);

    if (delta < this->delta_min()) {
      converged = true;
      this->exit_solver(it, ConvergenceReason::CONVERGED_TR_DELTA);
    }

    return converged;
  }

  virtual bool check_local_convergence(const SizeType &it,
                                       const SizeType &it_success,
                                       const SizeType &level,
                                       const Scalar &delta,
                                       const LocalSolveType &solve_type) {
    if (this->check_iter_convergence(it, it_success, level, solve_type)) {
      return true;
    } else if (delta < this->delta_min()) {
      return true;
    }

    return this->criticality_measure_termination(level);
  }

  virtual bool check_iter_convergence(const SizeType &it,
                                      const SizeType &it_success,
                                      const SizeType &level,
                                      const LocalSolveType &solve_type) {
    // coarse one
    if (level == 0 && (it_success >= this->max_sucessful_coarse_it() ||
                       it >= this->max_coarse_it())) {
      return true;
    }
    // every other level
    else if (level > 0 && solve_type == PRE_SMOOTHING) {
      if (it >= this->pre_smoothing_steps() ||
          it_success >= this->max_sucessful_smoothing_it()) {
        return true;
      }
    } else if (level > 0 && solve_type == POST_SMOOTHING) {
      if (it >= this->post_smoothing_steps() ||
          it_success >= this->max_sucessful_smoothing_it()) {
        return true;
      }
    }

    return false;
  }

  virtual Scalar criticality_measure(const SizeType &level) = 0;
  virtual bool recursion_termination_smoothness(const Vector &g_restricted,
                                                const Vector &g_coarse,
                                                const SizeType & /*level*/) = 0;

  virtual bool criticality_measure_termination(const SizeType &level) {
    // this should be ideally different norm of grad on previous iteration
    Scalar atol_level =
        (level == this->n_levels() - 1)
            ? this->atol()
            : std::min(this->atol(), this->grad_smoothess_termination() *
                                         this->memory_.gnorm[level + 1]);
    return (this->memory_.gnorm[level] < atol_level) ? true : false;
  }

  ///////////////////////////////////////////// TR-based stuff
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual bool solve_qp_subproblem(const SizeType &level, const bool &flg) = 0;

  virtual void initialize_local_solve(const SizeType & /*level*/,
                                      const LocalSolveType & /*solve_type*/) {}

  virtual bool delta_update(const Scalar &rho, const SizeType &level,
                            const Vector &s_global) = 0;

  virtual Scalar get_pred(const SizeType &level) = 0;

 protected:
  RMTRLevelMemory<Matrix, Vector> memory_;
  MultilevelDerivEval<Matrix, Vector, CONSISTENCY_LEVEL> ml_derivs_;

  bool init_;
};

}  // namespace utopia

#endif  // UTOPIA_RMTR_BASE_HPP
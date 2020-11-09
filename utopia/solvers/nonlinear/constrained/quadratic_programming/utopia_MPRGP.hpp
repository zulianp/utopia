#ifndef UTOPIA_CONSTRAINT_QP_MPRGP
#define UTOPIA_CONSTRAINT_QP_MPRGP

#include <string>
#include "utopia_Algorithms.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_For.hpp"
#include "utopia_QPSolver.hpp"

namespace utopia {

template <class Matrix, class Vector>
class MPGRP final : public OperatorBasedQPSolver<Matrix, Vector> {
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;
  using Solver = utopia::LinearSolver<Matrix, Vector>;
  using Super = utopia::OperatorBasedQPSolver<Matrix, Vector>;

 public:
  using Super::solve;
  using Super::update;

  MPGRP() : eps_eig_est_(1e-1), power_method_max_it_(10) {}

  void read(Input &in) override {
    OperatorBasedQPSolver<Matrix, Vector>::read(in);
    in.get("eig_comp_tol", eps_eig_est_);
    in.get("power_method_max_it", power_method_max_it_);
  }

  void print_usage(std::ostream &os) const override {
    OperatorBasedQPSolver<Matrix, Vector>::print_usage(os);
    this->print_param_usage(os, "eig_comp_tol", "double",
                            "Tolerance of eigen solver.", "1e-1");
    this->print_param_usage(
        os, "power_method_max_it", "int",
        "Maximum number of iterations used inside of power method.", "10");
  }

  MPGRP *clone() const override { return new MPGRP(*this); }

  void update(const Operator<Vector> &A) override {
    UTOPIA_TRACE_REGION_BEGIN("MPGRP::update");

    const auto layout_rhs = row_layout(A);
    if (!initialized_ || !layout_rhs.same(layout_)) {
      init_memory(layout_rhs);
    }

    UTOPIA_TRACE_REGION_END("MPGRP::update");
  }

  bool solve(const Operator<Vector> &A, const Vector &rhs,
             Vector &sol) override {
    UTOPIA_TRACE_REGION_BEGIN("MPRGP::solve(...)");

    if (this->has_empty_bounds()) {
      this->fill_empty_bounds(layout(rhs));
    } else {
      assert(this->get_box_constraints().valid(layout(rhs)));
    }

    auto &box = this->get_box_constraints();

    this->update(A);

    // as it is not clear ATM, how to apply preconditioner, we use it at least
    // to obtain initial guess
    if (precond_) {
      // this is unconstrained step
      precond_->apply(rhs, sol);
      // projection to feasible set
      this->make_iterate_feasible(sol);
    }

    bool ok = aux_solve(A, rhs, sol, box);

    UTOPIA_TRACE_REGION_END("MPRGP::solve(...)");
    return ok;
  }

  void set_eig_comp_tol(const Scalar &eps_eig_est) {
    eps_eig_est_ = eps_eig_est;
  }

 private:
  bool aux_solve(const Operator<Vector> &A, const Vector &rhs, Vector &x,
                 const BoxConstraints<Vector> &constraints) {
    // UTOPIA_NO_ALLOC_BEGIN("MPRGP");
    // //cudaProfilerStart();

    // Scalar r_norm0 = norm2(rhs);

    const auto &&ub = constraints.upper_bound();
    const auto &&lb = constraints.lower_bound();

    if (this->verbose()) {
      this->init_solver("MPGRP comm_size: " + std::to_string(rhs.comm().size()),
                        {"it", "|| g ||"});
    }

    const Scalar gamma = 1.0;
    const Scalar alpha_bar = 1.95 / this->get_normA(A);
    Scalar pAp, beta_beta, fi_fi, gp_dot;

    SizeType it = 0;
    bool converged = false;
    Scalar gnorm;

    Scalar alpha_cg, alpha_f, beta_sc;

    // this->get_projection(x, *lb, *ub, Ax);

    assert(lb);
    assert(ub);

    // utopia::out() <<x.comm().size() << " " << rhs.comm().size() << " " <<
    // lb->comm().size() << " "
    //           << ub->comm().size() << std::endl;

    this->project(*lb, *ub, x);
    // x = Ax;

    A.apply(x, Ax);
    g = Ax - rhs;

    // std::cout<<"loc_size-lb: "<< local_size(*lb).get(0) << "  \n";
    // std::cout<<"loc_size-ub: "<< local_size(*ub).get(0) << "  \n";
    // std::cout<<"loc_size-g: "<< local_size(g).get(0) << "  \n";
    // std::cout<<"loc_size-x: "<< local_size(x).get(0) << "  \n";
    // std::cout<<"loc_size-fi: "<< local_size(fi).get(0) << "  \n";

    this->get_fi(x, g, *lb, *ub, fi);

    // std::cout<<"----------------------- \n";
    // exit(0);

    this->get_beta(x, g, *lb, *ub, beta);

    gp = fi + beta;
    p = fi;

    dots(beta, beta, beta_beta, fi, fi, fi_fi);

    while (!converged) {
      if (beta_beta <= (gamma * gamma * fi_fi)) {
        A.apply(p, Ap);

        dots(p, Ap, pAp, g, p, gp_dot);

        alpha_cg = gp_dot / pAp;
        y = x - alpha_cg * p;
        alpha_f = get_alpha_f(x, p, *lb, *ub, help_f1, help_f2);

        if (alpha_cg <= alpha_f) {
          x = y;
          g = g - alpha_cg * Ap;
          this->get_fi(x, g, *lb, *ub, fi);
          beta_sc = dot(fi, Ap) / pAp;
          p = fi - beta_sc * p;
        } else {
          x = x - alpha_f * p;
          g = g - alpha_f * Ap;
          this->get_fi(x, g, *lb, *ub, fi);

          help_f1 = x - (alpha_bar * fi);
          this->project(help_f1, *lb, *ub, x);

          A.apply(x, Ax);
          g = Ax - rhs;
          this->get_fi(x, g, *lb, *ub, p);
        }
      } else {
        A.apply(beta, Abeta);
        alpha_cg = dot(g, beta) / dot(beta, Abeta);
        x = x - alpha_cg * beta;
        g = g - alpha_cg * Abeta;

        this->get_fi(x, g, *lb, *ub, p);
      }

      this->get_fi(x, g, *lb, *ub, fi);
      this->get_beta(x, g, *lb, *ub, beta);

      gp = fi + beta;

      dots(beta, beta, beta_beta, fi, fi, fi_fi, gp, gp, gnorm);

      gnorm = std::sqrt(gnorm);
      it++;

      if (this->verbose()) {
        PrintInfo::print_iter_status(it, {gnorm});
      }

      converged = this->check_convergence(it, gnorm, 1, 1);
      // converged = (it > this->max_it() || gnorm < std::min(0.1,
      // std::sqrt(r_norm0)) * r_norm0 ) ? true : false;
    }

    // //cudaProfilerStop();
    // UTOPIA_NO_ALLOC_END();
    return true;
  }

 public:
  void get_fi(const Vector &x, const Vector &g, const Vector &lb,
              const Vector &ub, Vector &fi) const {
    assert(!empty(fi));

    {
      auto d_lb = const_local_view_device(lb);
      auto d_ub = const_local_view_device(ub);
      auto d_x = const_local_view_device(x);
      auto d_g = const_local_view_device(g);
      auto d_fi = local_view_device(fi);

      parallel_for(local_range_device(fi), UTOPIA_LAMBDA(const SizeType i) {
        // read all
        const Scalar li = d_lb.get(i);
        const Scalar ui = d_ub.get(i);
        const Scalar xi = d_x.get(i);
        const Scalar gi = d_g.get(i);

        d_fi.set(i, (li < xi && xi < ui) ? gi : Scalar(0.0));
      });
    }

    // {

    //     auto d_lb = const_device_view(lb);
    //     auto d_ub = const_device_view(ub);
    //     auto d_x  = const_device_view(x);
    //     auto d_g  = const_device_view(g);

    //     parallel_each_write(fi, UTOPIA_LAMBDA(const SizeType i) -> Scalar
    //     {
    //         Scalar li = d_lb.get(i);
    //         Scalar ui = d_ub.get(i);
    //         Scalar xi = d_x.get(i);
    //         Scalar gi = d_g.get(i);

    //         if(li < xi && xi < ui){
    //             return gi;
    //         }
    //         else{
    //             return 0.0;
    //         }

    //     });
    // }
  }

  Scalar get_alpha_f(const Vector &x, const Vector &p, const Vector &lb,
                     const Vector &ub, Vector &help_f1, Vector &help_f2) const {
    assert(!empty(help_f1));
    assert(!empty(help_f2));

    {
      auto d_lb = const_local_view_device(lb);
      auto d_ub = const_local_view_device(ub);
      auto d_x = const_local_view_device(x);
      auto d_p = const_local_view_device(p);

      auto h1 = local_view_device(help_f1);
      auto h2 = local_view_device(help_f2);

      parallel_for(local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
        // read all for quantities
        const Scalar li = d_lb.get(i);
        const Scalar ui = d_ub.get(i);
        const Scalar xi = d_x.get(i);
        const Scalar pi = d_p.get(i);

        // write both helpers
        h1.set(i, (pi > 0) ? ((xi - li) / pi) : Scalar(1e15));
        h2.set(i, (pi < 0) ? ((xi - ui) / pi) : Scalar(1e15));
      });
    }

    return multi_min(help_f1, help_f2);

    // {
    //     auto d_lb = const_device_view(lb);
    //     auto d_ub = const_device_view(ub);
    //     auto d_x  = const_device_view(x);
    //     auto d_p  = const_device_view(p);

    //     parallel_each_write(help_f1, UTOPIA_LAMBDA(const SizeType i) ->
    //     Scalar
    //     {
    //         Scalar li = d_lb.get(i);
    //         Scalar xi = d_x.get(i);
    //         Scalar pi = d_p.get(i);

    //         if(pi > 0)
    //         {
    //             return (xi-li)/pi;
    //         }
    //         else
    //         {
    //             return 1e15;
    //         }
    //     });

    //     parallel_each_write(help_f2, UTOPIA_LAMBDA(const SizeType i) ->
    //     Scalar
    //     {
    //         Scalar ui = d_ub.get(i);
    //         Scalar xi = d_x.get(i);
    //         Scalar pi = d_p.get(i);

    //         if(pi < 0)
    //         {
    //             return (xi-ui)/pi;
    //         }
    //         else
    //         {
    //             return 1e15;
    //         }

    //     });
    // }

    // return multi_min(help_f1, help_f2);
  }

  void get_beta(const Vector &x, const Vector &g, const Vector &lb,
                const Vector &ub, Vector &beta) const {
    assert(!empty(beta));

    {
      auto d_lb = const_local_view_device(lb);
      auto d_ub = const_local_view_device(ub);
      auto d_x = const_local_view_device(x);
      auto d_g = const_local_view_device(g);
      auto d_beta = local_view_device(beta);

      parallel_for(local_range_device(beta), UTOPIA_LAMBDA(const SizeType i) {
        const Scalar li = d_lb.get(i);
        const Scalar ui = d_ub.get(i);
        const Scalar xi = d_x.get(i);
        const Scalar gi = d_g.get(i);

        const Scalar val =
            (device::abs(li - xi) < 1e-14)
                ? device::min(0.0, gi)
                : ((device::abs(ui - xi) < 1e-14) ? device::max(0.0, gi) : 0.0);

        d_beta.set(i, val);
      });
    }

    // {
    //     auto d_lb = const_device_view(lb);
    //     auto d_ub = const_device_view(ub);
    //     auto d_x  = const_device_view(x);
    //     auto d_g  = const_device_view(g);

    //     parallel_each_write(beta, UTOPIA_LAMBDA(const SizeType i) -> Scalar
    //     {
    //         Scalar li = d_lb.get(i);
    //         Scalar ui = d_ub.get(i);
    //         Scalar xi = d_x.get(i);
    //         Scalar gi = d_g.get(i);

    //         if(device::abs(li -  xi) < 1e-14)
    //         {
    //             return device::min(0.0, gi);
    //         }
    //         else if(device::abs(ui -  xi) < 1e-14)
    //         {
    //             return device::max(0.0, gi);
    //         }
    //         else
    //         {
    //             return 0.0;
    //         }
    //     });
    // }
  }

 private:
  Scalar get_normA(const Operator<Vector> &A) {
    // Super simple power method to estimate the biggest eigenvalue
    assert(!empty(help_f2));
    help_f2.set(1.0);

    SizeType it = 0;
    bool converged = false;
    Scalar gnorm, lambda = 0.0, lambda_old;

    while (!converged) {
      help_f1 = help_f2;
      A.apply(help_f1, help_f2);
      help_f2 = Scalar(1.0 / Scalar(norm2(help_f2))) * help_f2;

      lambda_old = lambda;

      A.apply(help_f2, help_f1);
      lambda = dot(help_f2, help_f1);

      fi = help_f2 - help_f1;
      gnorm = norm2(fi);

      converged = ((gnorm < eps_eig_est_) ||
                   (std::abs(lambda_old - lambda) < eps_eig_est_) ||
                   it > power_method_max_it_)
                      ? true
                      : false;

      it = it + 1;
    }

    if (this->verbose())
      utopia::out() << "Power method converged in " << it
                    << " iterations. Largest eig: " << lambda << "  \n";

    return lambda;
  }

 public:
  void init_memory(const Layout &layout) override {
    OperatorBasedQPSolver<Matrix, Vector>::init_memory(layout);

    fi.zeros(layout);
    beta.zeros(layout);
    gp.zeros(layout);
    p.zeros(layout);
    y.zeros(layout);
    Ap.zeros(layout);
    Abeta.zeros(layout);
    Ax.zeros(layout);
    g.zeros(layout);
    help_f1.zeros(layout);
    help_f2.zeros(layout);

    initialized_ = true;
    layout_ = layout;
  }

  void set_preconditioner(
      const std::shared_ptr<Preconditioner<Vector> > &precond) override {
    precond_ = precond;
  }

 private:
  Vector fi, beta, gp, p, y, Ap, Abeta, Ax, g, help_f1, help_f2;

  Scalar eps_eig_est_;
  SizeType power_method_max_it_;

  bool initialized_{false};
  Layout layout_;

  std::shared_ptr<Preconditioner<Vector> > precond_;
};
}  // namespace utopia

#endif  // UTOPIA_CONSTRAINT_QP_MPRGP

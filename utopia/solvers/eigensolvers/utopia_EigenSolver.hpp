#ifndef UTOPIA_EIGEN_SOLVER_HPP
#define UTOPIA_EIGEN_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Smoother.hpp"

namespace utopia {

template <typename Matrix, typename Vector>
class EigenSolver : public virtual Clonable, public Configurable {
 public:
  using Scalar = typename utopia::Traits<Vector>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;

  EigenSolver()
      : number_of_eigenvalues_(1),
        portion_of_spectrum_("smallest_real"),
        max_it_(40000),
        tol_(1e-12)

  {}

  ~EigenSolver() override = default;

  virtual void number_of_eigenvalues(const SizeType &number_of_eigenvalues) {
    number_of_eigenvalues_ = number_of_eigenvalues;
  }

  virtual const SizeType &number_of_eigenvalues() const {
    return number_of_eigenvalues_;
  }

  virtual void max_it(const SizeType &max_it) { max_it_ = max_it; }

  virtual const SizeType &max_it() const { return max_it_; }

  virtual void tol(const Scalar &tol) { tol_ = tol; }

  virtual const Scalar &tol() const { return tol_; }

  virtual const bool &verbose() const { return verbose_; }

  virtual void verbose(const bool &verbose) { verbose_ = verbose; }

  virtual const std::string &portion_of_spectrum() const {
    return portion_of_spectrum_;
  }

  void read(Input &in) override {
    in.get("max_it", max_it_);
    in.get("tol", tol_);
    in.get("verbose", verbose_);

    in.get("number_of_eigenvalues", number_of_eigenvalues_);
    in.get("portion_of_spectrum", portion_of_spectrum_);
  }

  void print_usage(std::ostream &os) const override {
    this->print_param_usage(os, "max_it", "int",
                            "Maximum number of iterations.", "1000");
    this->print_param_usage(os, "tol", "real", "Absolute tolerance.", "0.25");
    this->print_param_usage(os, "verbose", "bool", "Verbose flag.", "false");

    this->print_param_usage(os, "number_of_eigenvalues", "int",
                            "Number of eigenvalues to be computed.", "1");
    this->print_param_usage(os, "portion_of_spectrum", "string",
                            "Define portion of spectrum of interest.",
                            "smallest_real");
  }

  virtual void portion_of_spectrum(const std::string &type) = 0;

  virtual bool solve(const Matrix &A) = 0;
  virtual bool solve(const Matrix &A, const Matrix &B) = 0;

  virtual bool print_eigenpairs() = 0;
  virtual void get_eigenpairs(const SizeType &i, Scalar &iegr, Scalar &eigi,
                              Vector &vr, Vector &vi) = 0;
  virtual void get_real_eigenpair(const SizeType &i, Scalar &iegr,
                                  Vector &vr) = 0;

  EigenSolver *clone() const override = 0;

 private:
  SizeType number_of_eigenvalues_;
  std::string portion_of_spectrum_;

  SizeType max_it_;
  Scalar tol_;
  bool verbose_{false};
};

}  // namespace utopia

#endif  // UTOPIA_EIGEN_SOLVER_HPP
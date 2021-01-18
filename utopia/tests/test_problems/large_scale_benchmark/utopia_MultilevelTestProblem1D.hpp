#ifndef UTOPIA_MULTILEVEL_TEST_PROBLEM_1D_HPP
#define UTOPIA_MULTILEVEL_TEST_PROBLEM_1D_HPP

#include "utopia.hpp"

namespace utopia {
template <typename Matrix, typename Vector>
class MultilevelTestProblemBase {
 public:
  using SizeType = typename Traits<Vector>::SizeType;
  using Scalar = typename Traits<Vector>::SizeType;

  MultilevelTestProblemBase(const SizeType &n_levels = 2,
                            const SizeType &n_coarse = 10,
                            const bool remove_bc = false)
      : n_levels_(n_levels), n_coarse_(n_coarse), remove_bc_(remove_bc) {
    n_dofs_.resize(n_levels);
    transfers_.resize(n_levels - 1);
    level_functions_.resize(n_levels);
  }

  virtual ~MultilevelTestProblemBase() {
    level_functions_.clear();
    transfers_.clear();
    n_dofs_.clear();
  }

  SizeType n_levels() const { return n_levels_; }

  SizeType n_coarse() const { return n_coarse_; }

  SizeType n_dofs(const SizeType &level) const { return n_dofs_[level]; }

  void n_dofs(const SizeType &level, const SizeType &dofs_l) {
    n_dofs_[level] = dofs_l;
  }

  bool remove_bc() const { return remove_bc_; }

  const std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > &get_transfer()
      const {
    return transfers_;
  }

  const std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >
      &get_functions() const {
    return level_functions_;
  }

 private:
  SizeType n_levels_;
  SizeType n_coarse_;
  bool remove_bc_;
  std::vector<SizeType> n_dofs_;

 protected:
  std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > transfers_;
  std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >
      level_functions_;
};

template <class Matrix, class Vector, class ProblemType>
class MultiLevelTestProblem1D final
    : public MultilevelTestProblemBase<Matrix, Vector> {
 public:
  using Traits = utopia::Traits<Vector>;
  using Scalar = typename Traits::Scalar;
  using SizeType = typename Traits::SizeType;
  using Comm = typename Traits::Communicator;

  MultiLevelTestProblem1D(const SizeType &n_levels,
                          const SizeType &n_coarse_elements,
                          const bool remove_bc = false,
                          const bool scale_by_h = false)
      : MultilevelTestProblemBase<Matrix, Vector>(n_levels, n_coarse_elements,
                                                  remove_bc) {
    assert(n_coarse_elements > 0);
    assert(n_levels > 0);

    this->n_dofs(0, n_coarse_elements + 1);

    // std::cout<<"n_coarse_elements: "<< n_coarse_elements << "  \n";

    Comm comm;

    for (SizeType i = 1; i < n_levels; ++i) {
      this->n_dofs(i, (this->n_dofs(i - 1) - 1) * 2 + 1);
    }

    for (SizeType i = 0; i < n_levels - 1; ++i) {
      Scalar h = scale_by_h ? (1. / (this->n_dofs(i - 1) - 1)) : 1.0;

      const auto n_coarse = this->n_dofs(i);
      const auto n_fine = this->n_dofs(i + 1);

      // std::cout<<"n_fine: "<< n_fine << "  n_ccoarse: "<< n_coarse << "  \n";

      Matrix I;
      I.sparse(
          layout(comm, Traits::decide(), Traits::decide(), n_fine, n_coarse), 2,
          2);
      // std::cout<<"n_coarse: "<< n_coarse << "  n_fine: "<< n_fine << "  \n";

      {
        Write<Matrix> w_(I, utopia::GLOBAL_INSERT);

        auto r = row_range(I);

        SizeType j = std::floor(r.begin() / 2.);

        SizeType reminder = r.begin() % 2;
        SizeType r_begin = r.begin() + reminder;

        if (reminder) {
          if (j + 1 < n_coarse) {
            I.set(r.begin(), j, 0.5 / h);
            I.set(r.begin(), j + 1, 0.5 / h);
          }

          ++j;
        }

        for (auto k = r_begin; k < r.end(); k += 2, ++j) {
          I.set(k, j, 1. / h);

          if (j + 1 < n_coarse) {
            auto kp1 = k + 1;

            if (r.inside(kp1)) {
              I.set(kp1, j, 0.5 / h);
              I.set(kp1, j + 1, 0.5 / h);
            }
          }
        }
      }

      // disp(I, "I");

      if (this->remove_bc() && i == n_levels - 2) {
        Write<Matrix> w_(I);
        auto rr = row_range(I);

        if (rr.inside(0)) {
          I.set(0, 0, 0.);
        }

        auto last_node_h = size(I).get(0) - 1;
        auto last_node_H = size(I).get(1) - 1;
        if (rr.inside(last_node_h)) {
          I.set(last_node_h, last_node_H, 0.);
        }
      }

      // approx of projection ...
      auto h_fine = (1. / (this->n_dofs(i + 1) - 1));
      auto h_coarse = (1. / (this->n_dofs(i) - 1));
      Matrix P = 1. / h_coarse * h_fine * transpose(I);
      // disp(P, "P");

      {
        Write<Matrix> w_(P, utopia::GLOBAL_INSERT);
        auto rP = row_range(P);

        auto last_node_H = size(P).get(1) - 1.0;
        auto last_node_h = size(P).get(0) - 1;

        if (rP.inside(last_node_h)) {
          P.set(last_node_h, last_node_H - 1, 0.0);
          P.set(last_node_h, last_node_H, 1.0);
        }

        if (rP.inside(0)) {
          P.set(0, 1, 0.);
          P.set(0, 0, 1.);
        }
      }

      // this->transfers_[i] = std::make_shared<IPTransfer<Matrix, Vector> >(
      // std::make_shared<Matrix>(I)); this->transfers_[i] =
      // std::make_shared<IPRTransfer<Matrix, Vector> >(
      // std::make_shared<Matrix>(I));
      this->transfers_[i] = std::make_shared<IPRTransfer<Matrix, Vector> >(
          std::make_shared<Matrix>(I), std::make_shared<Matrix>(P));
    }

    for (auto l = 0; l < n_levels; l++) {
      auto fun = std::make_shared<ProblemType>(this->n_dofs(l));
      this->level_functions_[l] = fun;
    }
  }
};
}  // namespace utopia
#endif
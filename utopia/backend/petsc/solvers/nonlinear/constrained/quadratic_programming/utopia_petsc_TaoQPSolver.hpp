#ifndef UTOPIA_PETSC_TAO_QP_SOLVEWR_HPP
#define UTOPIA_PETSC_TAO_QP_SOLVEWR_HPP

#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_QPSolver.hpp"

#include "utopia_petsc_Factorization.hpp"

#include <memory>

namespace utopia {

    template <class Matrix, class Vector>
    class TaoQPSolver {};

    template <>
    class TaoQPSolver<PetscMatrix, PetscVector> final : public QPSolver<PetscMatrix, PetscVector> {
        using Scalar = utopia::Traits<PetscVector>::Scalar;
        using SizeType = utopia::Traits<PetscVector>::SizeType;
        typedef utopia::LinearSolver<PetscMatrix, PetscVector> LinearSolver;
        typedef utopia::QPSolver<PetscMatrix, PetscVector> Super;
        using BoxConstraints = utopia::BoxConstraints<PetscVector>;

    public:
        TaoQPSolver(const std::shared_ptr<LinearSolver> &linear_solver =
                        std::make_shared<Factorization<PetscMatrix, PetscVector>>());

        ~TaoQPSolver() override;
        TaoQPSolver *clone() const override;
        bool apply(const PetscVector &rhs, PetscVector &sol) override;
        void set_linear_solver(const std::shared_ptr<LinearSolver> &linear_solver);
        void tao_type(const std::string &type);
        void read(Input &in) override;

        Scalar atol() const override;
        Scalar rtol() const override;
        Scalar stol() const override;

        SizeType max_it() const override;
        bool verbose() const override;

        void atol(const Scalar &atol) override;
        void rtol(const Scalar &rtol) override;
        void stol(const Scalar &stol) override;
        void max_it(const SizeType &max_it) override;
        void verbose(const bool &verbose) override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_TAO_QP_SOLVEWR_HPP

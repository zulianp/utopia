#include "utopia_petsc_TaoQPSolver.hpp"

#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

#include "utopia_petsc_TaoQPSolver.hpp"
#include "utopia_petsc_TaoSolver.hpp"

#include "utopia_QuadraticFunction.hpp"
#include "utopia_make_unique.hpp"

#include "petsctao.h"

#include <cassert>

namespace utopia {

    class TaoQPSolver<PetscMatrix, PetscVector>::Impl {
    public:
        typedef utopia::LinearSolver<PetscMatrix, PetscVector> LinearSolver;

        Impl(const std::shared_ptr<LinearSolver> &linear_solver)
            : tao(utopia::make_unique<TaoSolver<PetscMatrix, PetscVector>>(linear_solver)) {
            tao->set_type(TAOGPCG);
        }

        std::unique_ptr<TaoSolver<PetscMatrix, PetscVector>> tao;
    };

    TaoQPSolver<PetscMatrix, PetscVector>::TaoQPSolver(const std::shared_ptr<LinearSolver> &linear_solver)
        : impl_(utopia::make_unique<Impl>(linear_solver)) {}

    TaoQPSolver<PetscMatrix, PetscVector>::~TaoQPSolver() = default;

    TaoQPSolver<PetscMatrix, PetscVector> *TaoQPSolver<PetscMatrix, PetscVector>::clone() const {
        m_utopia_warning_once("TaoQPSolver * TaoQPSolver::clone() not implemented properly");
        // FIXME
        // exception-safe cloning
        auto cloned_tao = std::unique_ptr<TaoSolver<PetscMatrix, PetscVector>>(impl_->tao->clone());
        auto cloned = utopia::make_unique<TaoQPSolver>();
        cloned->impl_->tao = std::move(cloned_tao);
        return cloned.release();
    }

    bool TaoQPSolver<PetscMatrix, PetscVector>::apply(const PetscVector &rhs, PetscVector &sol) {
        assert(this->has_operator());

        impl_->tao->atol(this->atol());
        impl_->tao->stol(this->stol());
        impl_->tao->rtol(this->rtol());
        impl_->tao->max_it(this->max_it());

        impl_->tao->set_box_constraints(this->get_box_constraints());

        QuadraticFunction<PetscMatrix, PetscVector> fun(std::make_shared<PetscMatrix>(*this->get_operator()),
                                                        std::make_shared<PetscVector>(rhs));

        return impl_->tao->solve(fun, sol);
    }

    void TaoQPSolver<PetscMatrix, PetscVector>::tao_type(const std::string &type) { impl_->tao->set_type(type); }

    void TaoQPSolver<PetscMatrix, PetscVector>::read(Input &in) { impl_->tao->read(in); }

    void TaoQPSolver<PetscMatrix, PetscVector>::set_linear_solver(const std::shared_ptr<LinearSolver> &linear_solver) {
        impl_->tao->set_linear_solver(linear_solver);
    }

    std::shared_ptr<utopia::LinearSolver<PetscMatrix, PetscVector>>
    TaoQPSolver<PetscMatrix, PetscVector>::linear_solver() const {
        return impl_->tao->linear_solver();
    }

    typename TaoQPSolver<PetscMatrix, PetscVector>::Scalar TaoQPSolver<PetscMatrix, PetscVector>::atol() const {
        return impl_->tao->atol();
    }

    typename TaoQPSolver<PetscMatrix, PetscVector>::Scalar TaoQPSolver<PetscMatrix, PetscVector>::rtol() const {
        return impl_->tao->rtol();
    }

    typename TaoQPSolver<PetscMatrix, PetscVector>::Scalar TaoQPSolver<PetscMatrix, PetscVector>::stol() const {
        return impl_->tao->stol();
    }

    typename TaoQPSolver<PetscMatrix, PetscVector>::SizeType TaoQPSolver<PetscMatrix, PetscVector>::max_it() const {
        return impl_->tao->max_it();
    }

    bool TaoQPSolver<PetscMatrix, PetscVector>::verbose() const { return impl_->tao->verbose(); }

    void TaoQPSolver<PetscMatrix, PetscVector>::atol(const Scalar &atol) { impl_->tao->atol(atol); }

    void TaoQPSolver<PetscMatrix, PetscVector>::rtol(const Scalar &rtol) { impl_->tao->rtol(rtol); }

    void TaoQPSolver<PetscMatrix, PetscVector>::stol(const Scalar &stol) { impl_->tao->stol(stol); }

    void TaoQPSolver<PetscMatrix, PetscVector>::max_it(const SizeType &max_it) { impl_->tao->max_it(max_it); }

    void TaoQPSolver<PetscMatrix, PetscVector>::verbose(const bool &verbose) { impl_->tao->verbose(verbose); }

}  // namespace utopia

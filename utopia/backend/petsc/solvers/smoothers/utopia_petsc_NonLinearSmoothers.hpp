#ifndef UTOPIA_PETSC_NONLINEAR_SMOOTHER_HPP
#define UTOPIA_PETSC_NONLINEAR_SMOOTHER_HPP

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Smoother.hpp"

#include "utopia_petsc_SNES.hpp"
#include "utopia_petsc_SNESFunction.hpp"

#include <petsc/private/snesimpl.h>
#include "petscsnes.h"

#include "utopia_NonlinearSolverInterfaces.hpp"

// TODO:: check nonlinear preconditioners ...
namespace utopia {

    template <class Matrix, class Vector>
    class NonLinearGaussSeidel<Matrix, Vector, PETSC> final : public SNESSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        typedef utopia::SNESSolver<Matrix, Vector> SNESSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

    public:
        NonLinearGaussSeidel(const std::shared_ptr<LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>())
            : SNESSolver(linear_solver) {
            this->set_snes_type("ngs");
        }

    private:
        void set_snes_options(SNES &snes,
                              const Scalar &atol = SNESSolver::atol(),
                              const Scalar &rtol = SNESSolver::rtol(),
                              const Scalar &stol = SNESSolver::stol(),
                              const SizeType &max_it = SNESSolver::max_it()) override {
            SNESSolver::set_snes_options(snes, atol, rtol, stol, max_it);

            // we need to allocate hessian for coloring computation

            PetscBool assembled;
            MatAssembled(snes->jacobian_pre, &assembled);

            if (!assembled) SNESComputeJacobian(snes, snes->vec_sol, snes->jacobian, snes->jacobian_pre);

            SNESLineSearch linesearch;
            SNESGetLineSearch(snes, &linesearch);
            SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC);
        }
    };

    template <class Matrix, class Vector>
    class NonLinearConjugateGradient<Matrix, Vector, PETSC> final : public SNESSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        typedef utopia::SNESSolver<Matrix, Vector> SNESSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

    public:
        NonLinearConjugateGradient(const std::shared_ptr<LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(),
                                   const std::vector<std::string> &update_types = {"FR", "PRP", "HS", "DY", "CD"})
            : SNESSolver(linear_solver), update_types(update_types) {
            this->set_snes_type("ncg");
            update_type_ = update_types.at(0);
        }

        void update_type(const std::string &update_type) {
            update_type_ = in_array(update_type, update_types) ? update_type : update_types.at(0);
        }

        std::string &update_type() { return update_type_; }

        void read(Input &in) override {
            SNESSolver::read(in);

            std::string update_type_aux;
            in.get("update_type", update_type_aux);

            // validation check
            this->update_type(update_type_aux);
        }

        void print_usage(std::ostream &os) const override {
            SNESSolver::print_usage(os);
            this->print_param_usage(os, "update_type", "string", "Choice of update type.", "FR");
        }

    private:
        void set_snes_options(SNES &snes,
                              const Scalar &atol = SNESSolver::atol(),
                              const Scalar &rtol = SNESSolver::rtol(),
                              const Scalar &stol = SNESSolver::stol(),
                              const SizeType &max_it = SNESSolver::max_it()) override {
            SNESSolver::set_snes_options(snes, atol, rtol, stol, max_it);

            if (update_type_ == "PRP")
                SNESNCGSetType(snes, SNES_NCG_PRP);
            else if (update_type_ == "HS")
                SNESNCGSetType(snes, SNES_NCG_HS);
            else if (update_type_ == "DY")
                SNESNCGSetType(snes, SNES_NCG_DY);
            else if (update_type_ == "CD")
                SNESNCGSetType(snes, SNES_NCG_CD);
            else
                SNESNCGSetType(snes, SNES_NCG_FR);

            // error oriented LS seems to work the best ...
            SNESLineSearch linesearch;
            SNESGetLineSearch(snes, &linesearch);
            SNESLineSearchSetType(linesearch, SNESLINESEARCHCP);
        }

    private:
        const std::vector<std::string> update_types;
        std::string update_type_;
    };

    // TODO:: put more options
    template <class Matrix, class Vector>
    class NonLinearGMRES<Matrix, Vector, PETSC> final : public SNESSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        typedef utopia::SNESSolver<Matrix, Vector> SNESSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

    public:
        NonLinearGMRES(const std::shared_ptr<LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>())
            : SNESSolver(linear_solver) {
            this->set_snes_type("ngmres");
        }

    private:
        void set_snes_options(SNES &snes,
                              const Scalar &atol = SNESSolver::atol(),
                              const Scalar &rtol = SNESSolver::rtol(),
                              const Scalar &stol = SNESSolver::stol(),
                              const SizeType &max_it = SNESSolver::max_it()) override {
            SNESSolver::set_snes_options(snes, atol, rtol, stol, max_it);

            // error oriented LS seems to work the best ...
            SNESLineSearch linesearch;
            SNESGetLineSearch(snes, &linesearch);
            SNESLineSearchSetType(linesearch, SNESLINESEARCHCP);
        }
    };

    // TODO:: put more options
    template <class Matrix, class Vector>
    class NonLinearAnderson<Matrix, Vector, PETSC> final : public SNESSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        typedef utopia::SNESSolver<Matrix, Vector> SNESSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

    public:
        NonLinearAnderson(const std::shared_ptr<LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>())
            : SNESSolver(linear_solver) {
            this->set_snes_type("anderson");
        }

    private:
        void set_snes_options(SNES &snes,
                              const Scalar &atol = SNESSolver::atol(),
                              const Scalar &rtol = SNESSolver::rtol(),
                              const Scalar &stol = SNESSolver::stol(),
                              const SizeType &max_it = SNESSolver::max_it()) override {
            SNESSolver::set_snes_options(snes, atol, rtol, stol, max_it);

            // error oriented LS seems to work the best ...
            SNESLineSearch linesearch;
            SNESGetLineSearch(snes, &linesearch);
            SNESLineSearchSetType(linesearch, SNESLINESEARCHCP);
        }
    };

    // TODO:: put more options
    template <class Matrix, class Vector>
    class NonLinearRichardson<Matrix, Vector, PETSC> final : public SNESSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        typedef utopia::SNESSolver<Matrix, Vector> SNESSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

    public:
        NonLinearRichardson(const std::shared_ptr<LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(),
                            const Scalar &alpha = 1.0)
            : SNESSolver(linear_solver), alpha_(alpha) {
            this->set_snes_type("nrichardson");
        }

        void dumping_parameter(const Scalar &alpha) { alpha_ = alpha; }

        Scalar dumping_parameter() { return alpha_; }

        void read(Input &in) override {
            SNESSolver::read(in);

            std::string SNES_type_aux_;
            in.get("dumping_parameter", alpha_);
        }

        void print_usage(std::ostream &os) const override {
            SNESSolver::print_usage(os);
            this->print_param_usage(
                os, "dumping_parameter", "real", "Dumping parameter used to dump step direction.", "1.0");
        }

    private:
        void set_snes_options(SNES &snes,
                              const Scalar &atol = SNESSolver::atol(),
                              const Scalar &rtol = SNESSolver::rtol(),
                              const Scalar &stol = SNESSolver::stol(),
                              const SizeType &max_it = SNESSolver::max_it()) override {
            SNESSolver::set_snes_options(snes, atol, rtol, stol, max_it);

            SNESLineSearch linesearch;
            SNESGetLineSearch(snes, &linesearch);
            SNESLineSearchSetType(linesearch, SNESLINESEARCHCP);

            // set damping
            SNESLineSearchSetDamping(linesearch, alpha_);
        }

    private:
        Scalar alpha_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_NONLINEAR_SMOOTHER_HPP

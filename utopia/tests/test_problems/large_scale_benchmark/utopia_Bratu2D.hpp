#ifndef UTOPIA_BRATU2D_HPP
#define UTOPIA_BRATU2D_HPP

#include "utopia.hpp"
#include "utopia_TestFunctions.hpp"
#include "utopia_petsc_Layout.hpp"

namespace utopia {
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class Bratu2D {};
}  // namespace utopia

#ifdef WITH_PETSC
#warning "this code should go in the petsc backend folder"

#include <petsc/private/snesimpl.h> /* For SNES_Solve event */
#include <petscdm.h>
#include <petscdmda.h>
#include <petscmatlab.h>
#include <petscsnes.h>

namespace utopia {
    using ParamsBratu2D = struct ParamsBratu2D;
    struct ParamsBratu2D {
        PetscReal lambda; /* test problem parameter */
        PetscErrorCode (*mms_solution)(ParamsBratu2D *, const DMDACoor2d *, PetscScalar *);
    };

    PetscErrorCode Bratu2DFormObjectiveLocal(DMDALocalInfo *, PetscScalar **, PetscReal *, ParamsBratu2D *);
    PetscErrorCode Bratu2DFormFunctionLocal(DMDALocalInfo *, PetscScalar **, PetscScalar **, ParamsBratu2D *);
    PetscErrorCode Bratu2DFormJacobianLocal(DMDALocalInfo *, PetscScalar **, Mat, Mat, ParamsBratu2D *);

    PetscErrorCode Bratu2DFormInitialGuess(DM, ParamsBratu2D *, Vec);
    PetscErrorCode Bratu2DMMSSolution(ParamsBratu2D *, const DMDACoor2d *, PetscScalar *);

    PetscErrorCode Bratu2DFormBCData(DM da, ParamsBratu2D *user, Vec BC_flag, Vec BC_value);

    template <typename Matrix, typename Vector>
    class Bratu2D<Matrix, Vector, PETSC> final : virtual public UnconstrainedExtendedTestFunction<Matrix, Vector>,
                                                 virtual public ConstrainedExtendedTestFunction<Matrix, Vector> {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        Bratu2D(const SizeType &n, const Scalar &lambda = 3.0) : n_(n), setup_(false) {
            application_context_.lambda = (lambda >= 0 && lambda <= 6.8) ? lambda : 3.4;
            application_context_.mms_solution = Bratu2DMMSSolution;

            this->create_DM();
            this->setup_SNES();
            this->setup_application_context();
            setup_ = true;

            auto vl = layout(snes_->vec_sol);

            this->constraints_ =
                make_box_constaints(std::make_shared<Vector>(vl, -9e9), std::make_shared<Vector>(vl, 0.45));

            // TD::fix this
            exact_sol_.zeros(vl);
        }

        Bratu2D(const DM &dm, const Scalar &lambda = 3.0) : setup_(false) {
            application_context_.lambda = (lambda >= 0 && lambda <= 6.8) ? lambda : 3.4;
            application_context_.mms_solution = Bratu2DMMSSolution;
            // application_context_.mms_forcing    = Bratu2DMMSForcing;

            da_ = dm;
            // necessary to provide reasonable global dimension
            DMDAGetInfo(da_,
                        nullptr,
                        &n_,
                        nullptr,
                        nullptr,
                        nullptr,
                        nullptr,
                        nullptr,
                        nullptr,
                        nullptr,
                        nullptr,
                        nullptr,
                        nullptr,
                        nullptr);

            this->setup_SNES();
            this->setup_application_context();
            setup_ = true;

            // TODO::find out exact solution - should be possible to compute
            // exact_sol_  = zeros(1, 0); //WHY initialize it at all?
        }

        ~Bratu2D() override {
            if (setup_) {
                DMDestroy(&da_);
                SNESDestroy(&snes_);
            }
        }

        bool gradient(const Vector &x, Vector &g) const override {
            // initialization of gradient vector...
            if (empty(g)) {
                g.zeros(layout(x));
                ;
            }

            SNESComputeFunction(snes_, raw_type(x), raw_type(g));

            // remove BC from gradient...

            return true;
        }

        bool hessian(const Vector &x, Matrix &hessian) const override {
            SNESComputeJacobian(snes_, raw_type(x), snes_->jacobian, snes_->jacobian);

            // yes, wrap would be nicer, but lets not use it ...
            convert(snes_->jacobian, hessian);

            // disp(hessian);
            // write("hessian.m", hessian);

            return true;
        }

        bool value(const Vector &x, typename Vector::Scalar &result) const override {
            SNESComputeObjective(snes_, raw_type(x), &result);
            return true;
        }

        void output_to_VTK(const Vector &x, const std::string file_name = "Bratu2D.vtk") {
            PetscViewer viewer;

            PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK);
            DMView(da_, viewer);
            PetscObjectSetName((PetscObject)raw_type(x), "x");
            VecView(raw_type(x), viewer);
            PetscViewerDestroy(&viewer);
        }

        Vector initial_guess() const override {
            Vector x_utopia;
            convert(snes_->vec_sol, x_utopia);
            return x_utopia;
        }

        const Vector &exact_sol() const override { return exact_sol_; }

        Scalar min_function_value() const override {
            // depends on the solution to which we converged to
            return -1.012;
        }

        std::string name() const override { return "Bratu2D"; }

        SizeType dim() const override { return n_ * n_; }

        bool exact_sol_known() const override { return true; }

        bool parallel() const override { return true; }

        void lambda(const Scalar &lambda) { application_context_.lambda = lambda; }

        Scalar lambda() const { return application_context_.lambda; }

    private:
        void create_DM() {
            DMDACreate2d(PETSC_COMM_WORLD,
                         DM_BOUNDARY_NONE,
                         DM_BOUNDARY_NONE,
                         DMDA_STENCIL_STAR,
                         n_,
                         n_,
                         PETSC_DECIDE,
                         PETSC_DECIDE,
                         1,
                         1,
                         nullptr,
                         nullptr,
                         &da_);
            DMSetUp(da_);
            DMDASetUniformCoordinates(da_, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        }

        void setup_SNES() {
            SNESCreate(PETSC_COMM_WORLD, &snes_);
            SNESSetFromOptions(snes_);
            SNESSetDM(snes_, da_);
            DMDASNESSetObjectiveLocal(da_, (DMDASNESObjective)Bratu2DFormObjectiveLocal, &application_context_);
            DMDASNESSetFunctionLocal(
                da_, INSERT_VALUES, (DMDASNESFunction)Bratu2DFormFunctionLocal, &application_context_);
            DMDASNESSetJacobianLocal(da_, (DMDASNESJacobian)Bratu2DFormJacobianLocal, &application_context_);

            // preallocate matrices/vectors
            DMCreateMatrix(da_, &snes_->jacobian);
            DMCreateMatrix(da_, &snes_->jacobian_pre);
            DMCreateGlobalVector(da_, &snes_->vec_sol);
            Bratu2DFormInitialGuess(da_, &application_context_, snes_->vec_sol);
        }

        void setup_application_context() {
            DMSetApplicationContext(da_, &application_context_);

            // PetscInt n_loc;
            // VecGetLocalSize(snes_->vec_sol, &n_loc);

            auto vl = layout(snes_->vec_sol);
            Vector bc_markers(vl, 0.0);
            Vector bc_values(vl, 0.0);

            Bratu2DFormBCData(da_, &application_context_, raw_type(bc_markers), raw_type(bc_values));
            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, bc_values);
        }

    private:
        SizeType n_{0};  // global size
        bool setup_;

        ParamsBratu2D application_context_{};
        DM da_{nullptr};
        SNES snes_{nullptr};

        Vector exact_sol_;
        PetscCommunicator comm_;
    };
}  // namespace utopia

#endif  // WITH_PETSC
#endif

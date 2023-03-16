#ifndef UTOPIA_NONL_ELLIPSE_2D_HPP
#define UTOPIA_NONL_ELLIPSE_2D_HPP

#include "utopia.hpp"
#include "utopia_TestFunctions.hpp"

namespace utopia {
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class NonEllipse2D {};
}  // namespace utopia

#ifdef UTOPIA_WITH_PETSC
#include <petsc/private/snesimpl.h> /* For SNES_Solve event */
#include <petscdm.h>
#include <petscdmda.h>
#include <petscmatlab.h>
#include <petscsnes.h>

namespace utopia {

    template <typename Matrix, typename Vector>
    class NonEllipse2D<Matrix, Vector, PETSC> final : virtual public UnconstrainedExtendedTestFunction<Matrix, Vector>,
                                                      virtual public ConstrainedExtendedTestFunction<Matrix, Vector> {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        NonEllipse2D(const SizeType &n, const SizeType &problem_type = 2)
            : n_(n), setup_(false), problem_type_(problem_type), lambda_(10.0), pi_(3.14159265358979323846) {
            this->create_DM();
            this->setup_SNES();

            if (problem_type_ == 1) {
                setup_problem1();
            } else if (problem_type_ == 2) {
                setup_problem2();
            } else {
                utopia_error("NonEllipse2D:: problem type not valid. \n");
            }

            setup_ = true;
        }

        NonEllipse2D(const DM &dm, const SizeType &problem_type = 2)
            : setup_(false), problem_type_(problem_type), lambda_(10.0), pi_(3.14159265358979323846) {
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
            // this->setup_application_context();

            if (problem_type_ == 1) {
                setup_problem1();
            } else if (problem_type_ == 2) {
                setup_problem2();
            } else {
                utopia_error("NonEllipse2D:: problem type not valid. \n");
            }

            setup_ = true;
        }

        ~NonEllipse2D() override {
            if (setup_) {
                DMDestroy(&da_);
                SNESDestroy(&snes_);
            }
        }

        void get_A_rhs(Matrix &A, Vector &rhs) const {
            A = A_no_bc_;
            convert(snes_->vec_rhs, rhs);
        }

        bool gradient(const Vector &x, Vector &g) const override {
            // initialization of gradient vector...
            if (empty(g)) {
                g.zeros(layout(x));
            }

            // Vector rhs;
            UTOPIA_NO_ALLOC_BEGIN("Ellipse::grad1");
            convert(snes_->vec_rhs, *A_help_);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("Ellipse::grad2");
            g = A_no_bc_ * x + *A_help_;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("Ellipse::grad3");
            *A_help_ = exp(x);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("Ellipse::grad4");
            // *A_help_ = lambda_*e_mul(x, *A_help_);
            *A_help_ = e_mul(x, *A_help_);
            *A_help_ *= lambda_;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("Ellipse::grad5");
            g += HxHy_ * *A_help_;
            UTOPIA_NO_ALLOC_END();

            remove_BC_contrib(g);

            return true;
        }

        bool hessian(const Vector &x, Matrix &hessian) const override {
            // YES, wrap is more effiicient, but we do not want to own matrix ....
            // as RMTR, needs to modify hessian ...
            // wrap(snes_->jacobian, hessian);
            // convert(snes_->jacobian, hessian);

            hessian = A_no_bc_;

            Vector x_p1 = x;
            x_p1.shift(1.0);
            // Vector exp_term = exp(x_p1);
            // Vector exp_term = x_p1*exp(x_p1);

            Vector exp_term = e_mul(x_p1, exp(x));

            exp_term *= HxHy_ * lambda_;

            hessian += Matrix(diag(exp_term));

            const std::vector<SizeType> &index = this->get_indices_related_to_BC();
            set_zero_rows(hessian, index, 1.);

            // disp(hessian);
            // exit(0);

            return true;
        }

        bool value(const Vector &x, typename Vector::Scalar &result) const override {
            Vector res1 = 0.0 * x;
            Vector res2;
            convert(snes_->vec_rhs, res2);

            // MatMult(snes_->jacobian, raw_type(x), raw_type(res1));
            MatMult(raw_type(A_no_bc_), raw_type(x), raw_type(res1));

            // Poisson term
            result = (0.5 * dot(res1, x)) + dot(res2, x);

            Vector exp_x = exp(x);
            // Vector exp_term = e_mul(x, exp_x) - exp_x;
            result = result + lambda_ * (HxHy_ * (dot(x, exp_x) - sum(exp_x)));

            return true;
        }

        void output_to_VTK(const Vector &x, const std::string file_name = "NonEllipse2D.vtk") {
            PetscViewer viewer;

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 14, 0)
            PetscViewerVTKOpen(PETSC_COMM_WORLD, file_name.c_str(), FILE_MODE_WRITE, &viewer);
#else
            PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK);
#endif
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

        std::string name() const override { return "NonlEllipse2D"; }

        SizeType dim() const override { return n_ * n_ * n_; }

        bool exact_sol_known() const override { return false; }

        bool parallel() const override { return true; }

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
            // DMDASetInterpolationType(da_, DMDA_Q0);
        }

        bool setup_SNES() {
            SNESCreate(PETSC_COMM_WORLD, &snes_);
            SNESSetFromOptions(snes_);
            SNESSetDM(snes_, da_);

            // preallocate matrices/vectors
            DMCreateMatrix(da_, &snes_->jacobian);
            DMCreateGlobalVector(da_, &snes_->vec_sol);
            DMCreateGlobalVector(da_, &snes_->vec_rhs);

            return false;
        }

    private:
        void setup_problem1() {
            PetscInt M, N;
            DMDAGetInfo(da_,
                        nullptr,
                        &M,
                        &N,
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
            PetscScalar Hx = 1.0 / (PetscReal)(M);
            PetscScalar Hy = 1.0 / (PetscReal)(N);
            HxHy_ = Hx * Hy;

            lambda_ = 10.0;

            this->build_rhs();
            this->build_init_guess();
            this->build_hessian();

            PetscInt n_local, n_global;
            VecGetLocalSize(snes_->vec_sol, &n_local);
            VecGetSize(snes_->vec_sol, &n_global);

            auto vl = layout(comm_, n_local, n_global);

            exact_sol_.zeros(vl);

            Vector bc_markers(vl, 0.0);
            Vector bc_values(vl, 0.0);

            this->form_BC_marker(bc_markers, bc_values);
            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, bc_values);

            const std::vector<SizeType> &index = this->get_indices_related_to_BC();

            Matrix Hessian;
            wrap(snes_->jacobian, Hessian);
            A_no_bc_ = Hessian;
            set_zero_rows(Hessian, index, 1.);

            A_help_ = make_unique<Vector>(vl, 0.0);

            // this->constraints_ = make_box_constaints(std::make_shared<Vector>(local_values(n_loc, -9e9)),
            //                                          std::make_shared<Vector>(local_values(n_loc, 0.55)));
        }

        void setup_problem2() {
            PetscInt M, N;
            DMDAGetInfo(da_,
                        nullptr,
                        &M,
                        &N,
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
            PetscScalar Hx = 1.0 / (PetscReal)(M);
            PetscScalar Hy = 1.0 / (PetscReal)(N);
            HxHy_ = Hx * Hy;

            lambda_ = 1.0;

            this->build_rhs();
            this->build_init_guess();
            this->build_hessian();

            // PetscInt n_loc;
            // VecGetLocalSize(snes_->vec_sol, &n_loc);

            auto vl = layout(snes_->vec_sol);

            exact_sol_.zeros(vl);

            Vector bc_markers(vl, 0.0);
            Vector bc_values(vl, 0.0);

            this->form_BC_marker(bc_markers, bc_values);
            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, bc_values);

            const std::vector<SizeType> &index = this->get_indices_related_to_BC();

            Matrix Hessian;
            wrap(snes_->jacobian, Hessian);
            A_no_bc_ = Hessian;
            set_zero_rows(Hessian, index, 1.);

            Vector lb(vl, 0.0);
            Vector ub(vl, 0.6);
            form_lb2(lb);

            A_help_ = make_unique<Vector>(vl, 0.0);

            // this->constraints_ = make_upper_bound_constraints(std::make_shared<Vector>(ub));
            this->constraints_ = make_box_constaints(std::make_shared<Vector>(lb), std::make_shared<Vector>(ub));
        }

        bool build_hessian() {
            PetscInt i, j, M, N, xm, ym, xs, ys, num, numi, numj;
            PetscScalar v[5], Hx, Hy, HydHx, HxdHy;
            MatStencil row, col[5];

            DMDAGetInfo(da_,
                        nullptr,
                        &M,
                        &N,
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
            Hx = 1.0 / (PetscReal)(M);
            Hy = 1.0 / (PetscReal)(N);
            HxdHy = Hx / Hy;
            HydHx = Hy / Hx;
            DMDAGetCorners(da_, &xs, &ys, nullptr, &xm, &ym, nullptr);

            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    row.i = i;
                    row.j = j;

                    if (i == 0 || j == 0 || i == M - 1 || j == N - 1) {
                        num = 0;
                        numi = 0;
                        numj = 0;
                        if (j != 0) {
                            v[num] = -HxdHy;
                            col[num].i = i;
                            col[num].j = j - 1;
                            num++;
                            numj++;
                        }
                        if (i != 0) {
                            v[num] = -HydHx;
                            col[num].i = i - 1;
                            col[num].j = j;
                            num++;
                            numi++;
                        }
                        if (i != M - 1) {
                            v[num] = -HydHx;
                            col[num].i = i + 1;
                            col[num].j = j;
                            num++;
                            numi++;
                        }
                        if (j != N - 1) {
                            v[num] = -HxdHy;
                            col[num].i = i;
                            col[num].j = j + 1;
                            num++;
                            numj++;
                        }
                        v[num] = ((PetscReal)(numj)*HxdHy + (PetscReal)(numi)*HydHx);
                        col[num].i = i;
                        col[num].j = j;
                        num++;
                        MatSetValuesStencil(snes_->jacobian, 1, &row, num, col, v, INSERT_VALUES);
                    } else {
                        v[0] = -HxdHy;
                        col[0].i = i;
                        col[0].j = j - 1;
                        v[1] = -HydHx;
                        col[1].i = i - 1;
                        col[1].j = j;
                        v[2] = 2.0 * (HxdHy + HydHx);
                        col[2].i = i;
                        col[2].j = j;
                        v[3] = -HydHx;
                        col[3].i = i + 1;
                        col[3].j = j;
                        v[4] = -HxdHy;
                        col[4].i = i;
                        col[4].j = j + 1;
                        MatSetValuesStencil(snes_->jacobian, 1, &row, 5, col, v, INSERT_VALUES);
                    }
                }
            }

            MatAssemblyBegin(snes_->jacobian, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(snes_->jacobian, MAT_FINAL_ASSEMBLY);

            return true;
        }

        void form_lb2(Vector &lb) {
            PetscInt i, j, mx, my, xm, ym, xs, ys;
            PetscScalar **array_marker;

            DMDAGetInfo(da_,
                        PETSC_IGNORE,
                        &mx,
                        &my,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE);
            DMDAGetCorners(da_, &xs, &ys, nullptr, &xm, &ym, nullptr);
            DMDAVecGetArray(da_, raw_type(lb), &array_marker);

            DM coordDA;
            Vec coordinates;
            DMDACoor2d **coords;

            DMGetCoordinateDM(da_, &coordDA);
            DMGetCoordinates(da_, &coordinates);
            DMDAVecGetArray(coordDA, coordinates, &coords);

            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    PetscScalar x = coords[j][i].x;
                    PetscScalar y = coords[j][i].y;

                    array_marker[j][i] =
                        -8.0 * (x - 7.0 / 16.0) * (x - 7.0 / 16.0) - 8.0 * (y - 7. / 16.) * (y - 7. / 16.) + 0.2;
                }
            }

            DMDAVecRestoreArrayDOF(da_, raw_type(lb), &array_marker);
            DMDAVecRestoreArray(coordDA, coordinates, &coords);
            VecAssemblyBegin(raw_type(lb));
            VecAssemblyEnd(raw_type(lb));
        }

        void form_BC_marker(Vector &bc_marker, Vector &bc_values) {
            PetscInt i, j, /*k,*/ mx, my, xm, ym, xs, ys;
            PetscScalar **array_marker;
            PetscScalar **array_values;

            DMDAGetInfo(da_,
                        PETSC_IGNORE,
                        &mx,
                        &my,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE);
            DMDAGetCorners(da_, &xs, &ys, nullptr, &xm, &ym, nullptr);
            DMDAVecGetArray(da_, raw_type(bc_marker), &array_marker);
            DMDAVecGetArray(da_, raw_type(bc_values), &array_values);

            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    if (i == 0 || j == 0 || i == mx - 1 || j == my - 1) {
                        array_marker[j][i] = 1.0;
                        if (problem_type_ == 1) {
                            array_values[j][i] = 0.0;
                        } else if (problem_type_ == 2) {
                            array_values[j][i] = 0.0;
                        }
                    } else {
                        if (problem_type_ == 1 || problem_type_ == 2) {
                            array_marker[j][i] = 0.0;
                            array_values[j][i] = 0.0;
                        }
                    }
                }
            }

            DMDAVecRestoreArrayDOF(da_, raw_type(bc_marker), &array_marker);
            VecAssemblyBegin(raw_type(bc_marker));
            VecAssemblyEnd(raw_type(bc_marker));

            DMDAVecRestoreArrayDOF(da_, raw_type(bc_values), &array_values);
            VecAssemblyBegin(raw_type(bc_values));
            VecAssemblyEnd(raw_type(bc_values));
        }

        void build_init_guess() {
            PetscInt i, j, /*k,*/ mx, my, xm, ym, xs, ys;
            PetscScalar **array;

            DMDAGetInfo(da_,
                        PETSC_IGNORE,
                        &mx,
                        &my,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE);
            DMDAGetCorners(da_, &xs, &ys, nullptr, &xm, &ym, nullptr);
            DMDAVecGetArray(da_, snes_->vec_sol, &array);

            DM coordDA;
            Vec coordinates;
            DMDACoor2d **coords;

            DMGetCoordinateDM(da_, &coordDA);
            DMGetCoordinates(da_, &coordinates);
            DMDAVecGetArray(coordDA, coordinates, &coords);

            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    if (i == 0.0 || j == 0.0 || i == mx - (1.0) || j == my - (1.0)) {
                        if (problem_type_ == 1) {
                            array[j][i] = 0.0;
                        } else if (problem_type_ == 2) {
                            array[j][i] = 0.0;
                        }
                    } else {
                        array[j][i] = 0.0;
                    }
                }
            }

            DMDAVecRestoreArrayDOF(da_, snes_->vec_sol, &array);
            DMDAVecRestoreArray(coordDA, coordinates, &coords);
            VecAssemblyBegin(snes_->vec_sol);
            VecAssemblyEnd(snes_->vec_sol);
        }

        void build_rhs() {
            PetscInt /*d, dof,*/ i, j, /*k,*/ mx, my, xm, ym, xs, ys;
            PetscScalar **array;
            PetscScalar Hx, Hy;

            DMDAGetInfo(da_,
                        PETSC_IGNORE,
                        &mx,
                        &my,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE);

            DMDAGetCorners(da_, &xs, &ys, nullptr, &xm, &ym, nullptr);
            Hx = 1.0 / (PetscReal)(mx);
            Hy = 1.0 / (PetscReal)(my);

            DM coordDA;
            Vec coordinates;
            DMDACoor2d **coords;

            DMGetCoordinateDM(da_, &coordDA);
            DMGetCoordinates(da_, &coordinates);
            DMDAVecGetArray(coordDA, coordinates, &coords);

            DMDAVecGetArray(da_, snes_->vec_rhs, &array);

            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    PetscScalar x = coords[j][i].x;
                    PetscScalar y = coords[j][i].y;

                    if (problem_type_ == 1) {
                        PetscScalar exp_term = ((x * x) - (x * x * x)) * std::sin(3.0 * pi_ * y);
                        PetscScalar term1 = ((9.0 * pi_ * pi_) + std::exp(exp_term)) * (x * x - x * x * x);
                        PetscScalar term2 = 6.0 * x - 2.0;
                        array[j][i] = -1.0 * (Hx * Hy) * ((term1 + term2) * std::sin(3.0 * pi_ * y));
                    } else if (problem_type_ == 2) {
                        PetscScalar exp_term = ((x * x) - (x * x * x)) * std::sin(3.0 * pi_ * y);
                        PetscScalar term1 = (9.0 * pi_ * pi_) + (std::exp(exp_term)) * (x * x - x * x * x);
                        PetscScalar term2 = 6.0 * x - 2.0;
                        array[j][i] = -1.0 * (Hx * Hy) * ((term1 + term2) * std::sin(3.0 * pi_ * x));
                    }
                }
            }

            DMDAVecRestoreArrayDOF(da_, snes_->vec_rhs, &array);
            DMDAVecRestoreArray(coordDA, coordinates, &coords);

            VecAssemblyBegin(snes_->vec_rhs);
            VecAssemblyEnd(snes_->vec_rhs);
        }

        void remove_BC_contrib(Vector &x) const {
            PetscInt i, j, /*k,*/ mx, my, xm, ym, xs, ys;
            PetscScalar **array;

            DMDAGetInfo(da_,
                        PETSC_IGNORE,
                        &mx,
                        &my,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE,
                        PETSC_IGNORE);
            DMDAGetCorners(da_, &xs, &ys, nullptr, &xm, &ym, nullptr);
            DMDAVecGetArray(da_, raw_type(x), &array);

            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    if (i == 0 || j == 0 || i == mx - 1 || j == my - 1) {
                        array[j][i] = 0.0;
                    }
                }
            }

            DMDAVecRestoreArrayDOF(da_, raw_type(x), &array);
            VecAssemblyBegin(raw_type(x));
            VecAssemblyEnd(raw_type(x));
        }

    private:
        SizeType n_;
        bool setup_;

        DM da_{nullptr};
        SNES snes_{nullptr};

        Vector exact_sol_;
        Matrix A_no_bc_;

        SizeType problem_type_;
        Scalar lambda_;
        Scalar pi_;

        Scalar HxHy_;

        std::unique_ptr<Vector> A_help_;
        PetscCommunicator comm_;
    };
}  // namespace utopia

#endif  // UTOPIA_WITH_PETSC
#endif  // UTOPIA_NONL_ELLIPSE_2D_HPP

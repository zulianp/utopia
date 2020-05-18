#include "utopia.hpp"
#include "utopia_TestFunctions.hpp"

namespace utopia {
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class Poisson3D {};
}  // namespace utopia

#ifdef WITH_PETSC
#include <petsc/private/snesimpl.h> /* For SNES_Solve event */
#include <petscdm.h>
#include <petscdmda.h>
#include <petscmatlab.h>
#include <petscsnes.h>

namespace utopia {

    template <typename Matrix, typename Vector>
    class Poisson3D<Matrix, Vector, PETSC> final : virtual public UnconstrainedExtendedTestFunction<Matrix, Vector>,
                                                   virtual public ConstrainedExtendedTestFunction<Matrix, Vector> {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        Poisson3D(const SizeType &n, const SizeType &problem_type = 1)
            : n_(n), setup_(false), problem_type_(problem_type) {
            if (problem_type_ > 1) {
                utopia_error("Poisson3D:: problem type not valid. \n");
            }

            this->create_DM();
            this->setup_SNES();
            this->setup_application_context();
            setup_ = true;

            PetscInt n_local, n_global;
            VecGetLocalSize(snes_->vec_sol, &n_local);
            VecGetSize(snes_->vec_sol, &n_global);
            auto vl = layout(comm_, n_local, n_global);

            this->constraints_ =
                make_box_constaints(std::make_shared<Vector>(vl, -9e9), std::make_shared<Vector>(vl, 0.45));
        }

        Poisson3D(const DM &dm) : setup_(false) {
            da_ = dm;
            // necessary to provide reasonable global dimension
            DMDAGetInfo(da_, 0, &n_, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

            this->setup_SNES();
            this->setup_application_context();
            setup_ = true;
        }

        ~Poisson3D() override {
            if (setup_) {
                DMDestroy(&da_);
                SNESDestroy(&snes_);
            }
        }

        bool get_rhs(Vector &rhs) const {
            convert(snes_->vec_rhs, rhs);
            return true;
        }

        void get_A_rhs(Matrix &A, Vector &rhs) const {
            A = A_no_bc_;
            convert(snes_->vec_rhs, rhs);
        }

        bool gradient(const Vector &x, Vector &g) const override {
            // initialization of gradient vector...
            if (empty(g)) {
                g.zeros(layout(x));
                ;
            }

            MatMultAdd(snes_->jacobian, raw_type(x), snes_->vec_rhs, raw_type(g));
            remove_BC_contrib(g);

            // disp(g);

            return true;
        }

        bool hessian(const Vector & /*x*/, Matrix &hessian) const override {
            // YES, wrap is more effiicient, but we do not want to own matrix ....
            // as RMTR, needs to modify hessian ...
            // wrap(snes_->jacobian, hessian);

            convert(snes_->jacobian, hessian);

            return true;
        }

        bool value(const Vector &x, typename Vector::Scalar &result) const override {
            Vector res1 = 0.0 * x;
            Vector res2;
            convert(snes_->vec_rhs, res2);

            // MatMult(snes_->jacobian, raw_type(x), raw_type(res1));
            MatMult(raw_type(A_no_bc_), raw_type(x), raw_type(res1));

            result = 0.5 * dot(res1, x) + dot(res2, x);

            return true;
        }

        void output_to_VTK(const Vector &x, const std::string file_name = "Poisson3D.vtk") {
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

        Scalar min_function_value() const override { return -1.013634375000014e+01; }

        std::string name() const override { return "Poisson3D"; }

        SizeType dim() const override { return n_ * n_ * n_; }

        bool exact_sol_known() const override { return true; }

        bool parallel() const override { return true; }

    private:
        void create_DM() {
            DMDACreate3d(PETSC_COMM_WORLD,
                         DM_BOUNDARY_NONE,
                         DM_BOUNDARY_NONE,
                         DM_BOUNDARY_NONE,
                         DMDA_STENCIL_STAR,
                         n_,
                         n_,
                         n_,
                         PETSC_DECIDE,
                         PETSC_DECIDE,
                         PETSC_DECIDE,
                         1,
                         1,
                         0,
                         0,
                         0,
                         &da_);
            DMSetUp(da_);
            DMDASetUniformCoordinates(da_, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            DMDASetInterpolationType(da_, DMDA_Q0);
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

        void setup_application_context() {
            this->build_rhs();
            this->build_init_guess();
            this->build_hessian();

            PetscInt n_local, n_global;
            VecGetLocalSize(snes_->vec_sol, &n_local);
            VecGetSize(snes_->vec_sol, &n_global);
            auto vl = layout(comm_, n_local, n_global);

            exact_sol_.values(vl, 0.0);
            this->build_exact_sol();

            Vector bc_markers(vl, 0.0);
            Vector bc_values(vl, 0.0);

            this->form_BC_marker(bc_markers, bc_values);
            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, bc_values);

            const std::vector<SizeType> &index = this->get_indices_related_to_BC();

            Matrix Hessian;
            wrap(snes_->jacobian, Hessian);

            A_no_bc_ = Hessian;

            set_zero_rows(Hessian, index, 1.);
        }

    private:
        bool build_hessian() {
            PetscErrorCode ierr;

            PetscInt dof, i, j, k, d, mx, my, mz, xm, ym, zm, xs, ys, zs, num, numi, numj, numk;
            PetscScalar v[7], Hx, Hy, Hz, HyHzdHx, HxHzdHy, HxHydHz;
            MatStencil row, col[7];

            ierr = DMDAGetInfo(da_, 0, &mx, &my, &mz, 0, 0, 0, &dof, 0, 0, 0, 0, 0);
            CHKERRQ(ierr);

            Hx = 1.0 / (PetscReal)(mx);
            Hy = 1.0 / (PetscReal)(my);
            Hz = 1.0 / (PetscReal)(mz);
            HyHzdHx = Hy * Hz / Hx;
            HxHzdHy = Hx * Hz / Hy;
            HxHydHz = Hx * Hy / Hz;

            ierr = DMDAGetCorners(da_, &xs, &ys, &zs, &xm, &ym, &zm);
            CHKERRQ(ierr);
            for (k = zs; k < zs + zm; k++) {
                for (j = ys; j < ys + ym; j++) {
                    for (i = xs; i < xs + xm; i++) {
                        for (d = 0; d < dof; d++) {
                            row.i = i;
                            row.j = j;
                            row.k = k;
                            row.c = d;
                            if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz - 1) {
                                num = 0;
                                numi = 0;
                                numj = 0;
                                numk = 0;
                                if (k != 0) {
                                    v[num] = -HxHydHz;
                                    col[num].i = i;
                                    col[num].j = j;
                                    col[num].k = k - 1;
                                    col[num].c = d;
                                    num++;
                                    numk++;
                                }
                                if (j != 0) {
                                    v[num] = -HxHzdHy;
                                    col[num].i = i;
                                    col[num].j = j - 1;
                                    col[num].k = k;
                                    col[num].c = d;
                                    num++;
                                    numj++;
                                }
                                if (i != 0) {
                                    v[num] = -HyHzdHx;
                                    col[num].i = i - 1;
                                    col[num].j = j;
                                    col[num].k = k;
                                    col[num].c = d;
                                    num++;
                                    numi++;
                                }
                                if (i != mx - 1) {
                                    v[num] = -HyHzdHx;
                                    col[num].i = i + 1;
                                    col[num].j = j;
                                    col[num].k = k;
                                    col[num].c = d;
                                    num++;
                                    numi++;
                                }
                                if (j != my - 1) {
                                    v[num] = -HxHzdHy;
                                    col[num].i = i;
                                    col[num].j = j + 1;
                                    col[num].k = k;
                                    col[num].c = d;
                                    num++;
                                    numj++;
                                }
                                if (k != mz - 1) {
                                    v[num] = -HxHydHz;
                                    col[num].i = i;
                                    col[num].j = j;
                                    col[num].k = k + 1;
                                    col[num].c = d;
                                    num++;
                                    numk++;
                                }
                                v[num] =
                                    (PetscReal)(numk)*HxHydHz + (PetscReal)(numj)*HxHzdHy + (PetscReal)(numi)*HyHzdHx;
                                col[num].i = i;
                                col[num].j = j;
                                col[num].k = k;
                                col[num].c = d;
                                num++;
                                ierr = MatSetValuesStencil(snes_->jacobian, 1, &row, num, col, v, INSERT_VALUES);
                                CHKERRQ(ierr);
                            } else {
                                v[0] = -HxHydHz;
                                col[0].i = i;
                                col[0].j = j;
                                col[0].k = k - 1;
                                col[0].c = d;
                                v[1] = -HxHzdHy;
                                col[1].i = i;
                                col[1].j = j - 1;
                                col[1].k = k;
                                col[1].c = d;
                                v[2] = -HyHzdHx;
                                col[2].i = i - 1;
                                col[2].j = j;
                                col[2].k = k;
                                col[2].c = d;
                                v[3] = 2.0 * (HyHzdHx + HxHzdHy + HxHydHz);
                                col[3].i = i;
                                col[3].j = j;
                                col[3].k = k;
                                col[3].c = d;
                                v[4] = -HyHzdHx;
                                col[4].i = i + 1;
                                col[4].j = j;
                                col[4].k = k;
                                col[4].c = d;
                                v[5] = -HxHzdHy;
                                col[5].i = i;
                                col[5].j = j + 1;
                                col[5].k = k;
                                col[5].c = d;
                                v[6] = -HxHydHz;
                                col[6].i = i;
                                col[6].j = j;
                                col[6].k = k + 1;
                                col[6].c = d;
                                ierr = MatSetValuesStencil(snes_->jacobian, 1, &row, 7, col, v, INSERT_VALUES);
                                CHKERRQ(ierr);
                            }
                        }
                    }
                }
            }
            ierr = MatAssemblyBegin(snes_->jacobian, MAT_FINAL_ASSEMBLY);
            CHKERRQ(ierr);
            ierr = MatAssemblyEnd(snes_->jacobian, MAT_FINAL_ASSEMBLY);
            CHKERRQ(ierr);

            return true;
        }

        void form_BC_marker(Vector &bc_marker, Vector &bc_values) {
            PetscInt d, dof, i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
            PetscScalar ****array_marker;
            PetscScalar ****array_values;

            DMDAGetInfo(da_, 0, &mx, &my, &mz, 0, 0, 0, &dof, 0, 0, 0, 0, 0);

            PetscScalar Hx = 1.0 / (PetscReal)(mx);
            PetscScalar Hy = 1.0 / (PetscReal)(my);
            PetscScalar Hz = 1.0 / (PetscReal)(mz);

            DMDAGetCorners(da_, &xs, &ys, &zs, &xm, &ym, &zm);
            DMDAVecGetArrayDOF(da_, raw_type(bc_marker), &array_marker);
            DMDAVecGetArrayDOF(da_, raw_type(bc_values), &array_values);

            for (k = zs; k < zs + zm; k++) {
                for (j = ys; j < ys + ym; j++) {
                    for (i = xs; i < xs + xm; i++) {
                        for (d = 0; d < dof; d++) {
                            if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz - 1) {
                                array_marker[k][j][i][d] = 1.0;

                                PetscScalar x = i * Hx;
                                PetscScalar y = j * Hy;
                                PetscScalar z = k * Hz;

                                if (problem_type_ == 0) {
                                    array_values[k][j][i][d] =
                                        (2. * x * (1. - x)) + (2. * y * (1. - y)) + (2. * z * (1. - z));
                                } else if (problem_type_ == 1) {
                                    array_values[k][j][i][d] = x * (1. - x) * y * (1. - y) * z * (1. - z);
                                }
                            } else {
                                array_marker[k][j][i][d] = 0.0;
                                array_values[k][j][i][d] = 0.0;
                            }
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
            PetscInt d, dof, i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
            PetscScalar ****array;

            DMDAGetInfo(da_, 0, &mx, &my, &mz, 0, 0, 0, &dof, 0, 0, 0, 0, 0);

            PetscScalar Hx = 1.0 / (PetscReal)(mx);
            PetscScalar Hy = 1.0 / (PetscReal)(my);
            PetscScalar Hz = 1.0 / (PetscReal)(mz);

            DMDAGetCorners(da_, &xs, &ys, &zs, &xm, &ym, &zm);
            DMDAVecGetArrayDOF(da_, snes_->vec_sol, &array);
            for (k = zs; k < zs + zm; k++) {
                for (j = ys; j < ys + ym; j++) {
                    for (i = xs; i < xs + xm; i++) {
                        for (d = 0; d < dof; d++) {
                            if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz - 1) {
                                PetscScalar x = i * Hx;
                                PetscScalar y = j * Hy;
                                PetscScalar z = k * Hz;

                                if (problem_type_ == 0) {
                                    array[k][j][i][d] = (2. * x * (1. - x)) + (2. * y * (1. - y)) + (2. * z * (1. - z));
                                } else if (problem_type_ == 1) {
                                    array[k][j][i][d] = x * (1. - x) * y * (1. - y) * z * (1. - z);
                                }
                            } else {
                                array[k][j][i][d] = 0.0;
                            }
                        }
                    }
                }
            }
            DMDAVecRestoreArrayDOF(da_, snes_->vec_sol, &array);
            VecAssemblyBegin(snes_->vec_sol);
            VecAssemblyEnd(snes_->vec_sol);
        }

        void build_rhs() {
            PetscInt d, dof, i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
            PetscScalar ****array;
            PetscScalar Hx, Hy, Hz;

            DMDAGetInfo(da_, 0, &mx, &my, &mz, 0, 0, 0, &dof, 0, 0, 0, 0, 0);

            DMDAGetCorners(da_, &xs, &ys, &zs, &xm, &ym, &zm);
            Hx = 1.0 / (PetscReal)(mx);
            Hy = 1.0 / (PetscReal)(my);
            Hz = 1.0 / (PetscReal)(mz);

            DMDAVecGetArrayDOF(da_, snes_->vec_rhs, &array);
            for (k = zs; k < zs + zm; k++) {
                for (j = ys; j < ys + ym; j++) {
                    for (i = xs; i < xs + xm; i++) {
                        for (d = 0; d < dof; d++) {
                            PetscScalar x1 = i * Hx;
                            PetscScalar x2 = j * Hy;
                            PetscScalar x3 = k * Hz;

                            if (problem_type_ == 0) {
                                array[k][j][i][d] = -12.0 * (Hx * Hy * Hz);
                            } else if (problem_type_ == 1) {
                                array[k][j][i][d] =
                                    (-(2. * x1 * x2 * (x1 - 1.) * (x2 - 1.)) - (2. * x1 * x3 * (x1 - 1.) * (x3 - 1.)) -
                                     (2. * x2 * x3 * (x2 - 1.) * (x3 - 1.))) *
                                    (Hx * Hy * Hz);
                            }
                        }
                    }
                }
            }

            DMDAVecRestoreArrayDOF(da_, snes_->vec_rhs, &array);
            VecAssemblyBegin(snes_->vec_rhs);
            VecAssemblyEnd(snes_->vec_rhs);
        }

        void build_exact_sol() {
            PetscInt d, dof, i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
            PetscScalar ****array;
            PetscScalar Hx, Hy, Hz;

            DMDAGetInfo(da_, 0, &mx, &my, &mz, 0, 0, 0, &dof, 0, 0, 0, 0, 0);

            DMDAGetCorners(da_, &xs, &ys, &zs, &xm, &ym, &zm);
            Hx = 1.0 / (PetscReal)(mx);
            Hy = 1.0 / (PetscReal)(my);
            Hz = 1.0 / (PetscReal)(mz);

            DMDAVecGetArrayDOF(da_, raw_type(exact_sol_), &array);
            for (k = zs; k < zs + zm; k++) {
                for (j = ys; j < ys + ym; j++) {
                    for (i = xs; i < xs + xm; i++) {
                        for (d = 0; d < dof; d++) {
                            PetscScalar x1 = i * Hx;
                            PetscScalar x2 = j * Hy;
                            PetscScalar x3 = k * Hz;

                            if (problem_type_ == 0) {
                                array[k][j][i][d] =
                                    (2. * x1 * (1. - x1)) + (2. * x2 * (1. - x2)) + (2. * x3 * (1. - x3));
                            } else if (problem_type_ == 1) {
                                array[k][j][i][d] = x1 * (1. - x1) * x2 * (1. - x2) * x3 * (1. - x3);
                            }
                        }
                    }
                }
            }

            DMDAVecRestoreArrayDOF(da_, raw_type(exact_sol_), &array);
            VecAssemblyBegin(raw_type(exact_sol_));
            VecAssemblyEnd(raw_type(exact_sol_));
        }

        void remove_BC_contrib(Vector &x) const {
            PetscInt d, dof, i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
            PetscScalar ****array;

            DMDAGetInfo(da_, 0, &mx, &my, &mz, 0, 0, 0, &dof, 0, 0, 0, 0, 0);
            DMDAGetCorners(da_, &xs, &ys, &zs, &xm, &ym, &zm);
            DMDAVecGetArrayDOF(da_, raw_type(x), &array);
            for (k = zs; k < zs + zm; k++) {
                for (j = ys; j < ys + ym; j++) {
                    for (i = xs; i < xs + xm; i++) {
                        for (d = 0; d < dof; d++) {
                            if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz - 1) {
                                array[k][j][i][d] = 0.0;
                            }
                        }
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
        PetscCommunicator comm_;
    };
}  // namespace utopia

#endif  // WITH_PETSC

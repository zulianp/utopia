#ifndef UTOPIA_MEMBRANE_2D_HPP
#define UTOPIA_MEMBRANE_2D_HPP

#include "utopia.hpp"
#include "utopia_TestFunctions.hpp"

namespace utopia {
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class Membrane2D {};
}  // namespace utopia

#ifdef UTOPIA_WITH_PETSC
#include <petsc/private/snesimpl.h> /* For SNES_Solve event */
#include <petscdm.h>
#include <petscdmda.h>
#include <petscmatlab.h>
#include <petscsnes.h>

namespace utopia {

    /**
     * @brief This problem was suggested as test problem for constrained QP solvers in Domoradova, Dostal, Projector
     * preconditioning for partially bound-constrained quadratic optimization and in Gratton, Sartaer, Toint, Numerical
     * Experience with a Recursive Trust-Region Method for Multilevel Nonlinear Optimization
     */
    template <typename Matrix, typename Vector>
    class Membrane2D<Matrix, Vector, PETSC> final : virtual public UnconstrainedExtendedTestFunction<Matrix, Vector>,
                                                    virtual public ConstrainedExtendedTestFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;
        using Layout = typename Traits::Layout;

        Membrane2D(const SizeType &n) : n_(n), setup_(false) {
            this->create_DM();
            this->setup_SNES();
            this->setup_problem();

            setup_ = true;
        }

        Membrane2D(const DM &dm) : setup_(false) {
            da_ = dm;
            // necessary to provide reasonable global dimension
            DMDAGetInfo(da_, 0, &n_, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

            this->setup_SNES();
            this->setup_problem();

            setup_ = true;
        }

        ~Membrane2D() override {
            if (setup_) {
                DMDestroy(&da_);
                SNESDestroy(&snes_);
            }
        }

        void get_A_rhs(Matrix &A) const { A = A_no_bc_; }

        bool gradient(const Vector &x, Vector &g) const override {
            // initialization of gradient vector...
            if (empty(g)) {
                g.zeros(layout(x));
                ;
            }

            g = A_no_bc_ * x + ones_;
            remove_BC_contrib(g);

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
            Vector res1 = A_no_bc_ * x;

            // this could be possibly optimized to get one reduction
            result = (0.5 * dot(res1, x)) + Scalar(HxHy_ * sum(x));
            return true;
        }

        void output_to_VTK(const Vector &x, const std::string file_name = "Membrane2D.vtk") {
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

        Scalar min_function_value() const override { return -1.490057781426608e-01; }

        std::string name() const override { return "Membrane2D"; }

        SizeType dim() const override { return n_ * n_; }

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

            return false;
        }

    private:
        void setup_problem() {
            this->build_init_guess();
            this->build_hessian();

            // PetscInt n_loc;
            // VecGetLocalSize(snes_->vec_sol, &n_loc);

            auto vl = layout(snes_->vec_sol);

            exact_sol_.zeros(vl);
            // exists, but not known
            // this->build_exact_sol();

            Vector bc_markers(vl, 0.0);
            Vector bc_values(vl, 0.0);

            this->build_ones(vl);

            this->form_BC_marker(bc_markers, bc_values);
            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, bc_values);

            const std::vector<SizeType> &index = this->get_indices_related_to_BC();

            Matrix Hessian;
            wrap(snes_->jacobian, Hessian);
            A_no_bc_ = Hessian;
            set_zero_rows(Hessian, index, 1.);

            Vector lb(vl, 0.0);
            this->form_lb(lb);

            this->constraints_ = make_lower_bound_constraints(std::make_shared<Vector>(lb));
        }

        void form_lb(Vector &lb) {
            PetscInt i, j, mx, my, xm, ym, xs, ys;
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
            DMDAVecGetArray(da_, raw_type(lb), &array);

            DM coordDA;
            Vec coordinates;
            DMDACoor2d **coords;

            DMGetCoordinateDM(da_, &coordDA);
            DMGetCoordinates(da_, &coordinates);
            DMDAVecGetArray(coordDA, coordinates, &coords);

            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    if (j == my - 1) {
                        // Constrained is computed as upper part of circle with the radius R=1, and the centre S=(1;
                        // 0.5; -1.3)
                        PetscScalar x = coords[j][i].x;
                        // PetscScalar y = coords[j][i].y;

                        PetscScalar c = ((x - 0.5) * (x - 0.5)) - 1.0 + (1.3 * 1.3);
                        PetscScalar b = 2.0 * 1.3;
                        array[j][i] = (-b + std::sqrt(b * b - 4.0 * c)) / 2.0;

                    } else {
                        array[j][i] = -9e9;
                    }
                }
            }

            DMDAVecRestoreArrayDOF(da_, raw_type(lb), &array);
            VecAssemblyBegin(raw_type(lb));
            VecAssemblyEnd(raw_type(lb));
        }

        void build_ones(const Layout &ls) {
            PetscInt mx, my;
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

            Hx = 1.0 / (PetscReal)(mx);
            Hy = 1.0 / (PetscReal)(my);
            HxHy_ = Hx * Hy;

            ones_.values(ls, Hx * Hy);
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

        void form_BC_marker(Vector &bc_marker, Vector &bc_values) {
            PetscInt i, j, mx, my, xm, ym, xs, ys;
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

            DM coordDA;
            Vec coordinates;
            DMDACoor2d **coords;

            DMGetCoordinateDM(da_, &coordDA);
            DMGetCoordinates(da_, &coordinates);
            DMDAVecGetArray(coordDA, coordinates, &coords);

            for (j = ys; j < ys + ym; j++) {
                for (i = xs; i < xs + xm; i++) {
                    // if (i==0 || j==0 || i==mx-1 || j==my-1)
                    if (j == 0) {
                        array_marker[j][i] = 1.0;

                        // PetscScalar x1 = coords[j][i].x;
                        // PetscScalar x2 = coords[j][i].y;

                        array_values[j][i] = 0.0;

                    } else {
                        array_marker[j][i] = 0.0;
                        array_values[j][i] = 0.0;
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
            PetscInt i, j, mx, my, xm, ym, xs, ys;
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
                    if (j == 0) {
                        array[j][i] = 0.0;
                    } else {
                        array[j][i] = 0.0;
                    }
                }
            }

            DMDAVecRestoreArrayDOF(da_, snes_->vec_sol, &array);
            VecAssemblyBegin(snes_->vec_sol);
            VecAssemblyEnd(snes_->vec_sol);
        }

        void remove_BC_contrib(Vector &x) const {
            PetscInt i, j, mx, my, xm, ym, xs, ys;
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
                    if (j == 0) {
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
        Vector ones_;
        Scalar HxHy_;
    };
}  // namespace utopia

#endif  // UTOPIA_WITH_PETSC
#endif  // UTOPIA_MEMBRANE_2D_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_TRILINOS
#ifdef UTOPIA_WITH_PETSC

// include edsl components
#include "utopia.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_petsc_trilinos.hpp"
#include "utopia_trilinos.hpp"

#include "utopia_petsc_CrsView.hpp"

#include "utopia_AlgebraicMultigrid.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_ILU.hpp"

#include "utopia_MPITimeStatistics.hpp"

#include "utopia_petsc_DILUAlgorithm.hpp"

#include "utopia_petsc_AdditiveCorrectionTransfer.hpp"

// std
#include <cmath>

namespace utopia {

    void convert_mm_matrix_to_petsc_bin(Input &in) {
        std::string path_in, path_out;
        in.get("in", path_in);
        in.get("out", path_out);

        TpetraMatrix mat;
        read(path_in, mat);

        PetscMatrix petsc_mat;
        convert(mat, petsc_mat);
        write(path_out, petsc_mat);
    }

    UTOPIA_REGISTER_APP(convert_mm_matrix_to_petsc_bin);

    void convert_mm_vector_to_petsc_bin(Input &in) {
        std::string path_in, path_out;
        in.get("in", path_in);
        in.get("out", path_out);

        TpetraVector mat;
        read(path_in, mat);

        PetscVector petsc_mat;
        convert(mat, petsc_mat);
        write(path_out, petsc_mat);
    }

    UTOPIA_REGISTER_APP(convert_mm_vector_to_petsc_bin);

    template <class Matrix, class Vector>
    class LSSolveApp : public Configurable {
    public:
        using Scalar = typename Traits<Vector>::Scalar;

        void read(Input &in) override {
            Matrix A;
            Vector x, b, oracle;

            MPITimeStatistics stats(x.comm());

            //////////////////////////////////////////
            stats.start();

            std::string path_A, path_b, path_oracle, path_output = "out.mm";
            bool use_amg = true;
            int block_size = 1;
            bool write_matlab = false;
            bool rescale_with_diag = false;
            bool use_ksp = false;
            bool use_ksp_smoother = true;
            bool convert_to_block_matrix = false;
            bool amg_as_preconditioner = false;

            in.require("A", path_A);
            in.require("b", path_b);
            in.get("oracle", path_oracle);
            in.get("out", path_output);
            in.get("use_amg", use_amg);
            in.get("block_size", block_size);
            in.get("write_matlab", write_matlab);
            in.get("rescale_with_diag", rescale_with_diag);
            in.get("use_ksp", use_ksp);
            in.get("use_ksp_smoother", use_ksp_smoother);
            in.get("convert_to_block_matrix", convert_to_block_matrix);
            in.get("amg_as_preconditioner", amg_as_preconditioner);

            if (use_ksp) {
                use_amg = false;
            }

            if (path_A.empty()) {
                utopia::err() << "[Error] A undefined!!!\n";
                return;
            }

            utopia::read(path_A, A);

            if (path_b.empty()) {
                utopia::err() << "[Warning] b is undefined, using vector of ones!!!\n";
                b.values(row_layout(A), 1.0);
            } else {
                utopia::read(path_b, b);
            }

            if (!path_oracle.empty()) {
                utopia::read(path_oracle, oracle);
            }

            stats.stop_collect_and_restart("read_files");

            if (rescale_with_diag) {
                Vector d = 1. / diag(A);
                A.diag_scale_left(d);
                b = e_mul(d, b);
            }

            if (convert_to_block_matrix) {
                Matrix A_temp;
                A.convert_to_mat_baij(block_size, A_temp);
                A = std::move(A_temp);
            }

            // utopia::out() << "Matrix::type() " << A.type() << '\n';

            if (write_matlab) {
                rename("a", A);
                write("A.m", A);

                rename("b", b);
                write("B.m", b);
            }

            x.zeros(row_layout(A));

            std::shared_ptr<LinearSolver<Matrix, Vector>> solver;

            if (use_amg) {
                std::shared_ptr<IterativeSolver<Matrix, Vector>> smoother;

                if (use_ksp_smoother) {
                    auto ksp = std::make_shared<KSPSolver<Matrix, Vector>>();
                    ksp->pc_type("ilu");
                    ksp->ksp_type("richardson");
                    smoother = ksp;
                } else {
                    smoother = std::make_shared<ILU<Matrix, Vector>>();
                }

                // auto smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
                std::shared_ptr<LinearSolver<Matrix, Vector>> coarse_solver;
                if (convert_to_block_matrix) {
                    auto gmres = std::make_shared<GMRES<Matrix, Vector>>("ilu");
                    gmres->factor_set_pivot_in_blocks(true);
                    coarse_solver = gmres;
                } else {
                    coarse_solver = std::make_shared<Factorization<Matrix, Vector>>();
                }

                // std::shared_ptr<MatrixAgglomerator<Matrix>> agglomerator;

                // if (block_size == 2) {
                //     agglomerator = std::make_shared<BlockAgglomerate<Matrix, 2>>();
                // } else if (block_size == 3) {
                //     agglomerator = std::make_shared<BlockAgglomerate<Matrix, 3>>();
                // } else if (block_size == 4) {
                //     agglomerator = std::make_shared<BlockAgglomerate<Matrix, 4>>();
                // } else {
                //     agglomerator = std::make_shared<Agglomerate<Matrix>>();
                // }

                auto amg = std::make_shared<AlgebraicMultigrid<Matrix, Vector>>(smoother, coarse_solver);

                InputParameters inner_params;
                inner_params.set("block_size", block_size);
                inner_params.set("use_simd", true);
                inner_params.set("use_line_search", true);
                coarse_solver->read(inner_params);
                smoother->read(inner_params);

                if (amg_as_preconditioner) {
                    amg->max_it(1);
                    auto krylov_solver = std::make_shared<BiCGStab<Matrix, Vector>>();
                    krylov_solver->set_preconditioner(amg);
                    solver = krylov_solver;
                } else {
                    solver = amg;
                }

            } else {
                if (use_ksp) {
                    auto ksp = std::make_shared<KSPSolver<Matrix, Vector>>();
                    ksp->factor_set_pivot_in_blocks(convert_to_block_matrix);
                    solver = ksp;
                } else {
                    // auto ksp = std::make_shared<KSPSolver<Matrix, Vector>>();
                    // auto ilu = std::make_shared<ILU<Matrix, Vector>>();

                    // InputParameters inner_params;
                    // inner_params.set("block_size", block_size);

                    // ilu->read(inner_params);
                    // ilu->max_it(1);

                    // ksp->set_preconditioner(ilu);
                    // solver = ksp;

                    solver = std::make_shared<ILU<Matrix, Vector>>();
                }
            }

            solver->read(in);

            stats.stop_collect_and_restart("read_settings");

            Vector r = b - A * x;
            Scalar r_norm = norm2(r);

            if (!mpi_world_rank()) {
                utopia::out() << "norm_residual (pre): " << r_norm << "\n";
                utopia::out() << "ndofs " << x.size() << std::endl;
            }

            solver->solve(A, b, x);

            stats.stop_collect_and_restart("solve");

            r = b - A * x;
            r_norm = norm2(r);

            if (!mpi_world_rank()) {
                utopia::out() << "norm_residual (post): " << r_norm << "\n";
            }

            write(path_output, x);

            stats.stop_collect_and_restart("write");
            stats.describe(utopia::out().stream());
        }
    };

    void ls_solve(Input &in) {
#ifdef UTOPIA_WITH_PETSC
        std::string backend = "petsc";
#else
        std::string backend = "trilinos";
#endif

        in.get("backend", backend);

        if (backend == "petsc") {
            LSSolveApp<PetscMatrix, PetscVector> app;
            app.read(in);
        } else {
            // AMG does not support Trilinos yet
            // LSSolveApp<TpetraMatrix, TpetraVector> app;
            // app.read(in);
        }
    }

    UTOPIA_REGISTER_APP(ls_solve);

}  // namespace utopia

#endif  // UTOPIA_WITH_TRILINOS
#endif  // UTOPIA_WITH_PETSC

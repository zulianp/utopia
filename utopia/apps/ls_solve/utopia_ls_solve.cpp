#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_TRILINOS
#ifdef UTOPIA_WITH_PETSC

// include edsl components
#include "utopia.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_petsc_trilinos.hpp"
#include "utopia_trilinos.hpp"

#include "utopia_MPITimeStatistics.hpp"

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

            in.get("A", path_A);
            in.get("b", path_b);
            in.get("oracle", path_oracle);
            in.get("out", path_output);

            if (path_A.empty()) {
                utopia::err() << "[Error] A undefined!!!\n";
                return;
            }

            if (path_b.empty()) {
                utopia::err() << "[Warning] b is undefined, using vector of ones!!!\n";
                b.values(row_layout(A), 1.0);
            } else {
                utopia::read(path_b, b);
            }

            utopia::read(path_A, A);

            if (!path_oracle.empty()) {
                utopia::read(path_oracle, oracle);
            }

            stats.stop_collect_and_restart("read_files");

            x.zeros(row_layout(A));

            KSPSolver<Matrix, Vector> solver;
            solver.read(in);

            // Factorization<Matrix, Vector> solver;

            stats.stop_collect_and_restart("read_settings");

            Vector r = b - A * x;
            Scalar r_norm = norm2(r);
            utopia::out() << "norm_residual (pre): " << r_norm << "\n";

            utopia::out() << "ndofs " << x.size() << std::endl;
            solver.solve(A, b, x);

            stats.stop_collect_and_restart("solve");

            r = b - A * x;
            r_norm = norm2(r);

            utopia::out() << "norm_residual (post): " << r_norm << "\n";

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
            LSSolveApp<TpetraMatrix, TpetraVector> app;
            app.read(in);
        }
    }

    UTOPIA_REGISTER_APP(ls_solve);

}  // namespace utopia

#endif  // UTOPIA_WITH_TRILINOS
#endif  // UTOPIA_WITH_PETSC
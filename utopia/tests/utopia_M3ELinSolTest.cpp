
#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_M3ELINSOL
#include "utopia.hpp"
#include "utopia_InputParameters.hpp"
#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

namespace utopia {

#ifdef UTOPIA_WITH_PETSC
    void amg_UTOPIA_WITH_PETSC() {
        PetscVector rhs, x;
        PetscMatrix A;

        const bool binwrite = false;
        const std::string data_path = Utopia::instance().get("data_path");
        const std::string folder = data_path + "/mg";
        const std::string sysfile = "system.txt";
        read(folder + "/rhs.bin", rhs);
        read(folder + "/A.bin", A);

        // const std::string folder = data_path + "/laplace/matrices_for_petsc";
        // read(folder + "/f_rhs", rhs);
        // read(folder + "/f_A", A);

        x = local_zeros(local_size(rhs));

        ASPAMG<PetscMatrix, PetscVector> amg;
        amg.solve(A, rhs, x);
        amg.print_system(binwrite, sysfile);  // Example on how to print a linear system to file in M3E's format

        double res_norm = norm2(rhs - A * x);
        utopia_test_assert(res_norm < 1e-8);
    }

#endif  // UTOPIA_WITH_PETSC

    // #ifdef UTOPIA_WITH_BLAS
    //     void amg_UTOPIA_WITH_BLAS()
    //     {
    //         BlasVectord rhs, x;
    //         BlasMatrixd A;

    //         const bool binwrite = false;
    //         const std::string data_path = Utopia::instance().get("data_path");
    //         const std::string folder = data_path + "/mg_blas";
    //         const std::string sysfile = "system.txt";

    //         read(folder + "/rhs.txt", rhs);
    //         read(folder + "/lhs.txt", A);

    //         x = local_zeros(local_size(rhs));

    //         ASPAMG<BlasMatrixd, BlasVectord> amg;
    //         if(!amg.import("ASPAMG", data_path + "/xml/default.xml")) {
    //             InputParameters in;
    //             in.set("TspMaxit", 200);
    //             amg.read(in);
    //         }

    //         amg.solve(A, rhs, x);
    //         amg.print_system(binwrite, sysfile); // Example on how to print a linear system to file in M3E's format

    //         double res_norm = norm2(rhs - A * x);
    //         utopia_test_assert(res_norm < 1e-8);

    //     }
    // #endif //UTOPIA_WITH_BLAS

    static void m3e_lin_sol() {
#ifdef UTOPIA_WITH_PETSC
        UTOPIA_RUN_TEST(amg_UTOPIA_WITH_PETSC);
#endif  // UTOPIA_WITH_PETSC

        // #ifdef UTOPIA_WITH_BLAS
        //         UTOPIA_RUN_TEST(amg_UTOPIA_WITH_BLAS);
        // #endif //UTOPIA_WITH_BLAS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(m3e_lin_sol);
}  // namespace utopia

#endif

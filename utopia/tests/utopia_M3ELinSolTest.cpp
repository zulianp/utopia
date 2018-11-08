#include "utopia_M3ELinSolTest.hpp"
#include "utopia_Base.hpp"
#include "utopia.hpp"
#include "utopia_InputParameters.hpp"
#include "utopia_ui.hpp"

namespace utopia {
#ifndef WITH_M3ELINSOL
	void run_m3e_lin_sol_test() {}
#else

#ifdef WITH_PETSC
	void amg_with_petsc()
	{
		DVectord  rhs, x;
		DSMatrixd A;
		
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

		ASPAMG<DSMatrixd, DVectord> amg;
		amg.solve(A, rhs, x);
		amg.print_system(binwrite, sysfile); // Example on how to print a linear system to file in M3E's format

		double res_norm = norm2(rhs - A * x);
		utopia_test_assert(res_norm < 1e-8);

	}

#endif //WITH_PETSC


#ifdef WITH_BLAS
	void amg_with_blas()
	{
		Vectord rhs, x;
		CRSMatrixd A;
		
		const bool binwrite = false;
		const std::string data_path = Utopia::instance().get("data_path");
		const std::string folder = data_path + "/mg_blas";
		const std::string sysfile = "system.txt";
		
		read(folder + "/rhs.txt", rhs);
		read(folder + "/lhs.txt", A);

		x = local_zeros(local_size(rhs));
	
		ASPAMG<CRSMatrixd, Vectord> amg;
		if(!amg.import("ASPAMG", data_path + "/json/default.json")) {
			InputParameters in;
			in.set("TspMaxit", 200);
			amg.read(in);
		}

		amg.solve(A, rhs, x);
		amg.print_system(binwrite, sysfile); // Example on how to print a linear system to file in M3E's format

		double res_norm = norm2(rhs - A * x);
		utopia_test_assert(res_norm < 1e-8);

	}
#endif //WITH_BLAS

	void run_m3e_lin_sol_test()
	{
		UTOPIA_UNIT_TEST_BEGIN("M3ELinSolTest");

#ifdef WITH_PETSC
		UTOPIA_RUN_TEST(amg_with_petsc);
#endif //WITH_PETSC

#ifdef WITH_BLAS
		UTOPIA_RUN_TEST(amg_with_blas);
#endif //WITH_BLAS

		UTOPIA_UNIT_TEST_END("M3ELinSolTest");
	}

#endif
}

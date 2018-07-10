#include "utopia_M3ELinSolTest.hpp"
#include "utopia_Base.hpp"
#include "utopia.hpp"

namespace utopia {
#ifndef WITH_M3ELINSOL
	void run_m3e_lin_sol_test() {}
#else

#ifdef WITH_PETSC
	void amg_with_petsc()
	{
		DVectord  rhs, x;
		DSMatrixd A;
		
		const std::string data_path = Utopia::instance().get("data_path");
		// const std::string folder = data_path + "/mg";
		// read(folder + "/rhs.bin", rhs);
		// read(folder + "/A.bin", A);

		const std::string folder = data_path + "/laplace/matrices_for_petsc";
		
		read(folder + "/f_rhs", rhs);
		read(folder + "/f_A", A);

		x = local_zeros(local_size(rhs));

		ASPAMG<DSMatrixd, DVectord> amg;
		amg.update(make_ref(A));
		amg.apply(rhs, x);

		double res_norm = norm2(rhs - A * x);
		utopia_test_assert(res_norm < 1e-8);

	}

#endif //WITH_PETSC

	void run_m3e_lin_sol_test()
	{
		UTOPIA_UNIT_TEST_BEGIN("M3ELinSolTest");

#ifdef WITH_PETSC
		UTOPIA_RUN_TEST(amg_with_petsc);
#endif //WITH_PETSC

		UTOPIA_UNIT_TEST_END("M3ELinSolTest");
	}

#endif
}
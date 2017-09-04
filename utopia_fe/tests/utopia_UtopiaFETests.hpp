#ifndef UTOPIA_FE_TESTS_HPP
#define UTOPIA_FE_TESTS_HPP

namespace libMesh {
	class LibMeshInit;
}

namespace utopia {
	void run_all_utopia_fe_tests(libMesh::LibMeshInit &init);
}

#endif //UTOPIA_FE_TESTS_HPP

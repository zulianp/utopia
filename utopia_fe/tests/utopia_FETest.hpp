#ifndef UTOPIA_FE_TEST_HPP
#define UTOPIA_FE_TEST_HPP

#include <string>

namespace libMesh {
	class LibMeshInit;
}

namespace utopia {
	void run_all_tests(libMesh::LibMeshInit &init);
	void run_test(const std::string &test_name, libMesh::LibMeshInit &init);
}

#endif //UTOPIA_FE_TEST_HPP
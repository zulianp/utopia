#ifndef UTOPIA_NON_LINEAR_ELASTICITY_TEST_HPP
#define UTOPIA_NON_LINEAR_ELASTICITY_TEST_HPP 

namespace libMesh {
	class LibMeshInit;
}

namespace utopia {
	void run_non_linear_elasticity_test(libMesh::LibMeshInit &init);
}

#endif //UTOPIA_NON_LINEAR_ELASTICITY_TEST_HPP

#ifndef UTOPIA_FE_ENERGY_ASSEMBLY_TEST_HPP
#define UTOPIA_FE_ENERGY_ASSEMBLY_TEST_HPP

namespace libMesh {
	class LibMeshInit;
}

namespace utopia 
{
	void run_energy_test(libMesh::LibMeshInit &init);
}

#endif //UTOPIA_FE_ENERGY_ASSEMBLY_TEST_HPP

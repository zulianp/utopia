#ifndef UTOPIA_INTERSECT_TEST_HPP
#define UTOPIA_INTERSECT_TEST_HPP

#include "utopia_fe_base.hpp"



namespace libMesh {
	class LibMeshInit;
}


namespace utopia {
	void run_intersect_test(libMesh::LibMeshInit &init);
}



#endif //UTOPIA_INTERSECT_TEST_HPP

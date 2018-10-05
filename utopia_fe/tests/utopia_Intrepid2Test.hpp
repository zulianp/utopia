#ifndef UTOPIA_INTREPID2_TEST_HPP
#define UTOPIA_INTREPID2_TEST_HPP

#include "utopia_fe_base.hpp"



namespace libMesh {
	class LibMeshInit;
}


namespace utopia {
	void run_intrepid2_test(libMesh::LibMeshInit &init);
}



#endif //UTOPIA_INTREPID2_TEST_HPP

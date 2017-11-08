#ifndef FEDSL_LEAST_SQUARES_EXAMPLES_HPP
#define FEDSL_LEAST_SQUARES_EXAMPLES_HPP

#include "utopia_fe_core.hpp"
#include "utopia_libmesh_Types.hpp"

namespace libMesh {
class LibMeshInit;
}

namespace utopia {

	LMDenseMatrix make_stress_strain_rel_tensor(const int dim, const double mu, const double lambda);
	void run_least_squares_examples(libMesh::LibMeshInit &init);
}

#endif //FEDSL_LEAST_SQUARES_EXAMPLES_HPP

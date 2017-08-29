#ifndef UTOPIA_SEMIGEOMETRIC_MULTIGRID_TEST_HPP
#define UTOPIA_SEMIGEOMETRIC_MULTIGRID_TEST_HPP 


#include "utopia_fe.hpp"
#include "utopia_assemble_volume_transfer.hpp"

// namespace libMesh {
// 	class LibMeshInit;
// }

namespace utopia {

	class MGTestProblem {
	public:
		typedef utopia::FESpace<LibMeshTraits<libMesh::Real> > FESpaceT;
		typedef std::shared_ptr<FESpaceT> FESpacePtr;
		typedef utopia::LibMeshFEContext<libMesh::LinearImplicitSystem> FEContextT;
		typedef std::shared_ptr<FEContextT> FEContextPtr;

		std::shared_ptr<libMesh::Mesh> coarse_mesh, fine_mesh;
		FESpacePtr coarse_space, fine_space;
		FEContextPtr coarse_context, fine_context;
		std::vector<std::shared_ptr<DSMatrixd>> interpolation_operators;

		DSMatrixd A;
		DVectord rhs;
	};

	void init_mg_test_problem(libMesh::LibMeshInit &init, MGTestProblem &problem);
	void run_semigeometric_multigrid_test(libMesh::LibMeshInit &init);
}

#endif //UTOPIA_SEMIGEOMETRIC_MULTIGRID_TEST_HPP

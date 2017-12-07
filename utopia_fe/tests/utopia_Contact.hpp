#ifndef UTOPIA_FE_CONTACT_HPP
#define UTOPIA_FE_CONTACT_HPP 

#include "utopia.hpp"
#include <vector>
#include <utility>

namespace libMesh {
	class MeshBase;
	class DofMap;

	namespace Parallel {
		class Communicator;
	}
}

namespace utopia {
	class Contact {
	public:
		bool init(const libMesh::Parallel::Communicator &lm_comm,
				  const std::shared_ptr<libMesh::MeshBase> &mesh,
				  const std::shared_ptr<libMesh::DofMap> &dof_map,
				  const double search_radius,
				  const std::vector<std::pair<int, int> > &contact_pair_tags,
				  unsigned int variable_number = 0);

		DVectord gap;
		DVectord weighted_gap;
		DVectord normals;
		DVectord inv_mass_vector;
		DVectord is_contact_node;

		DSMatrixd coupling;
		DSMatrixd inv_mass_matrix;
		DSMatrixd transfer_operator;
		DSMatrixd orthogonal_trafo;

		DSMatrixd complete_transformation;
	};
}

#endif //UTOPIA_FE_CONTACT_HPP


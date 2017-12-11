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
	class ContactParams {
	public:
		ContactParams()
		: search_radius(0.1), variable_number(0)
		{}

		double search_radius;
		std::vector<std::pair<int, int> > contact_pair_tags;
		unsigned int variable_number;
	};

	class Contact {
	public:
		Contact()
		: initialized(false)
		{}

		inline bool init(
			const std::shared_ptr<libMesh::MeshBase> &mesh,
			const std::shared_ptr<libMesh::DofMap> &dof_map,
			const ContactParams &params)
		{
			return init(mesh, dof_map, params.search_radius, params.contact_pair_tags, params.variable_number);
		}

		bool init(const std::shared_ptr<libMesh::MeshBase> &mesh,
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
		bool initialized;
	};
}

#endif //UTOPIA_FE_CONTACT_HPP


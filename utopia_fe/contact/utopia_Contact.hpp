#ifndef UTOPIA_FE_CONTACT_HPP
#define UTOPIA_FE_CONTACT_HPP 

#include "utopia.hpp"
#include "utopia_fe_base.hpp"

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
		: search_radius(0.1), variable_number(0), use_biorthogonal_basis(true)
		{}

		double search_radius;
		std::vector<std::pair<int, int> > contact_pair_tags;
		unsigned int variable_number;
		bool use_biorthogonal_basis;

		void describe(std::ostream &os) const;
	};

	class Contact {
	public:
		Contact()
		: initialized(false), has_contact_(false)
		{}

		bool init_no_contact(
			const std::shared_ptr<libMesh::MeshBase> &mesh,
			const std::shared_ptr<libMesh::DofMap> &dof_map);

		inline bool init(
			const std::shared_ptr<libMesh::MeshBase> &mesh,
			const std::shared_ptr<libMesh::DofMap> &dof_map,
			const ContactParams &params)
		{
			return init(mesh, dof_map, params.search_radius, params.contact_pair_tags, params.variable_number, params.use_biorthogonal_basis);
		}

		bool init(const std::shared_ptr<libMesh::MeshBase> &mesh,
				  const std::shared_ptr<libMesh::DofMap> &dof_map,
				  const double search_radius,
				  const std::vector<std::pair<int, int> > &contact_pair_tags,
				  unsigned int variable_number = 0,
				  const bool use_biorthogonal_basis = true);


		inline bool has_contact() const
		{
			return has_contact_;
		}

		UVector gap;
		UVector weighted_gap;
		UVector normals;
		UVector inv_mass_vector;
		UVector is_contact_node;

		USMatrix coupling;
		USMatrix inv_mass_matrix;
		USMatrix transfer_operator;
		USMatrix orthogonal_trafo;

		USMatrix complete_transformation;
		bool initialized;
		bool has_contact_;


		void print_debug_info();
	};
}

#endif //UTOPIA_FE_CONTACT_HPP


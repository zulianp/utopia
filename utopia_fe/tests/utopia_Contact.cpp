#include "utopia_Contact.hpp"

#include "utopia_assemble_contact.hpp"
#include "libmesh/parallel.h"
#include "libmesh/mesh_base.h"
#include "moonolith_communicator.hpp"
#include "moonolith_synched_describable.hpp"

#include <sstream>

namespace utopia {
	bool Contact::init(
		const std::shared_ptr<libMesh::MeshBase> &mesh,
		const std::shared_ptr<libMesh::DofMap> &dof_map,
		const double search_radius,
		const std::vector<std::pair<int, int> > &contact_pair_tags,
		unsigned int variable_number,
		const bool use_biorthogonal_basis)
	{

		moonolith::Communicator comm(mesh->comm().get());

		if(!assemble_contact(
			comm,
			mesh, 
			dof_map, 
			variable_number, 
			coupling, 
			orthogonal_trafo, 
			weighted_gap, 
			normals,
			is_contact_node, 
			search_radius,
			contact_pair_tags,
			use_biorthogonal_basis,
			false))
		{

		 	//something failed
			return false;
		}

		if(use_biorthogonal_basis) {

			DVectord d = sum(coupling, 1);

			inv_mass_vector = local_zeros(local_size(d));

			{
				Write<DVectord> w_(inv_mass_vector);

				each_read(d, [this](const SizeType i, const double value) {
					if(value < -1e-8) {
						std::cerr << "negative el for " << i << std::endl;
					}

					if(std::abs(value) > 1e-15) {
						this->inv_mass_vector.set(i, 1./value);
					} else {
						this->inv_mass_vector.set(i, 1.);
					}
				});
			}


			inv_mass_matrix = diag(inv_mass_vector);
		} else {
			assert(false && "implement me");

			// inv_mass_matrix = inv(boundary_mass_matrix);
		}


		transfer_operator = inv_mass_matrix * coupling;

		transfer_operator += local_identity(local_size(transfer_operator));
		gap = inv_mass_matrix * weighted_gap;
		complete_transformation = transfer_operator * orthogonal_trafo;

		initialized = true;
		return true;
	}

	void Contact::print_debug_info()
	{
		const double sum_T  = sum(transfer_operator);
		const double norm_g  = norm2(gap);
		const double norm_B  = norm2(coupling);
		const double norm_O  = norm2(orthogonal_trafo);
		const double norm_im = norm2(inv_mass_vector);

		std::stringstream ss;
		ss << "sum_T:   " << sum_T   << "\n";
		ss << "norm_g:  " << norm_g  << "\n";
		ss << "norm_B:  " << norm_B  << "\n";
		ss << "norm_O:  " << norm_O  << "\n";
		ss << "norm_im: " << norm_im << "\n";

		// static bool is_first = true;

		// if(is_first) {
		// 	DVectord t = sum(transfer_operator, 1);
		// 	write("t.m", t);
		// 	is_first = false;
		// }

		moonolith::Communicator comm;
		moonolith::root_describe(ss.str(), comm, std::cout);
	}
}
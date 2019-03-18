#include "utopia_Contact.hpp"

#include "utopia_assemble_contact.hpp"
#include "libmesh/parallel.h"
#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"
#include "moonolith_communicator.hpp"
#include "moonolith_synched_describable.hpp"

#include <sstream>

namespace utopia {

	void ContactParams::describe(std::ostream &os) const
	{
		os << "search_radius: " << search_radius << "\n";
		os << "variable_number: " << variable_number << "\n";
		os << "use_biorthogonal_basis: " << use_biorthogonal_basis << "\n";
		os << "master, slave:\n";
		for(const auto &p : contact_pair_tags) {
			os << p.first << ", " << p.second << "\n";
		}

		os << std::endl;
	}

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

			UVector d = sum(coupling, 1);

			inv_mass_vector = local_zeros(local_size(d));

			{
				Write<UVector> w_(inv_mass_vector);

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

		// write("B.m", coupling);
		// write("D_inv.m", inv_mass_matrix);
		// exit(0);


		transfer_operator = inv_mass_matrix * coupling;

		transfer_operator += local_identity(local_size(transfer_operator));
		gap = inv_mass_matrix * weighted_gap;
		complete_transformation = transfer_operator * orthogonal_trafo;

		has_contact_ = utopia::max(is_contact_node);

		initialized = true;
		return true;
	}

	bool Contact::init_no_contact(
		const std::shared_ptr<libMesh::MeshBase> &mesh,
		const std::shared_ptr<libMesh::DofMap> &dof_map)
	{
		auto n_local_dofs = dof_map->n_local_dofs();

		gap = local_values(n_local_dofs, 100000000);
		weighted_gap = gap;
		normals = local_zeros(n_local_dofs);

		inv_mass_vector = local_values(n_local_dofs, 1.);
		is_contact_node = local_zeros(n_local_dofs);

		coupling = local_identity(n_local_dofs, n_local_dofs);
		inv_mass_matrix = local_identity(n_local_dofs, n_local_dofs);
		transfer_operator = local_identity(n_local_dofs, n_local_dofs);
		orthogonal_trafo  = local_identity(n_local_dofs, n_local_dofs);
		complete_transformation = local_identity(n_local_dofs, n_local_dofs);
		initialized = true;
		has_contact_ = false;
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
		// 	UVector t = sum(transfer_operator, 1);
		// 	write("t.m", t);
		// 	is_first = false;
		// }

		moonolith::Communicator comm;
		moonolith::root_describe(ss.str(), comm, std::cout);
	}
}
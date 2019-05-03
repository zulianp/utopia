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
            coupling_,
            orthogonal_trafo_,
            weighted_gap_,
            normals_,
            is_contact_node_,
            search_radius,
            contact_pair_tags,
            use_biorthogonal_basis,
            false))
        {

             //something failed
            return false;
        }

        if(use_biorthogonal_basis) {

            UVector d = sum(coupling_, 1);

            inv_mass_vector_ = local_zeros(local_size(d));

            {
                Write<UVector> w_(inv_mass_vector_);

                each_read(d, [this](const SizeType i, const double value) {
                    if(value < -1e-8) {
                        std::cerr << "negative el for " << i << std::endl;
                    }

                    if(std::abs(value) > 1e-15) {
                        this->inv_mass_vector_.set(i, 1./value);
                    } else {
                        this->inv_mass_vector_.set(i, 1.);
                    }
                });
            }


            inv_mass_matrix_ = diag(inv_mass_vector_);
        } else {
            assert(false && "implement me");

            // inv_mass_matrix = inv(boundary_mass_matrix);
        }

        // write("B.m", coupling);
        // write("D_inv.m", inv_mass_matrix);
        // exit(0);


        transfer_operator_ = inv_mass_matrix_ * coupling_;

        transfer_operator_ += local_identity(local_size(transfer_operator_));
        gap_ = inv_mass_matrix_ * weighted_gap_;
        complete_transformation_ = transfer_operator_ * orthogonal_trafo_;

        has_contact_ = utopia::max(is_contact_node_);

        initialized_ = true;
        return true;
    }

    bool Contact::init_no_contact(
        const std::shared_ptr<libMesh::MeshBase> &mesh,
        const std::shared_ptr<libMesh::DofMap> &dof_map)
    {
        auto n_local_dofs = dof_map->n_local_dofs();

        gap_ = local_values(n_local_dofs, 100000000);
        weighted_gap_ = gap_;
        normals_ = local_zeros(n_local_dofs);

        inv_mass_vector_ = local_values(n_local_dofs, 1.);
        is_contact_node_ = local_zeros(n_local_dofs);

        coupling_ = local_identity(n_local_dofs, n_local_dofs);
        inv_mass_matrix_ = local_identity(n_local_dofs, n_local_dofs);
        transfer_operator_ = local_identity(n_local_dofs, n_local_dofs);
        orthogonal_trafo_  = local_identity(n_local_dofs, n_local_dofs);
        complete_transformation_ = local_identity(n_local_dofs, n_local_dofs);
        initialized_ = true;
        has_contact_ = false;
        return true;
    }

    void Contact::print_debug_info()
    {
        const double sum_T  = sum(transfer_operator_);
        const double norm_g  = norm2(gap_);
        const double norm_B  = norm2(coupling_);
        const double norm_O  = norm2(orthogonal_trafo_);
        const double norm_im = norm2(inv_mass_vector_);

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

    void Contact::couple(const UVector &in, UVector &out) const
    {
        out = transpose(complete_transformation_) * in;
    }
    
    void Contact::uncouple(const UVector &in, UVector &out) const
    {
        out = complete_transformation_ * in;
    }
    
    void Contact::couple(const USparseMatrix &in, USparseMatrix &out) const
    {
        out = transpose(complete_transformation_) * in * complete_transformation_;
    }
}
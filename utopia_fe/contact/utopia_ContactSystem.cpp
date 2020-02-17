#include "utopia_ContactSystem.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/numeric_vector.h"

namespace utopia {
    ContactSystem::ContactSystem(const std::shared_ptr<libMesh::EquationSystems> &equation_systems, const int main_system_num)
    : equation_systems_(equation_systems), main_system_num_(main_system_num), wear_coefficient_(1e-5)
    {
        init();
    }

    ContactSystem::~ContactSystem()
    {

    }

    void ContactSystem::init()
    {
        const int dim  = equation_systems_->get_mesh().mesh_dimension();
        auto &aux = equation_systems_->add_system<libMesh::LinearImplicitSystem>("aux");
        system_num_ = aux.number();

        auto &dof_map_main = equation_systems_->get_system(main_system_num_).get_dof_map();
        auto order = dof_map_main.variable_order(0);

        var_num_.push_back( aux.add_variable("inc_x", order, libMesh::LAGRANGE) );
        var_num_.push_back( aux.add_variable("inc_y", order, libMesh::LAGRANGE) );


        if(dim > 2)
            var_num_.push_back( aux.add_variable("inc_z", order, libMesh::LAGRANGE) );

        var_num_.push_back( aux.add_variable("vel_x", order, libMesh::LAGRANGE) );
        var_num_.push_back( aux.add_variable("vel_y", order, libMesh::LAGRANGE) );

        if(dim > 2)
            var_num_.push_back( aux.add_variable("vel_z", order, libMesh::LAGRANGE) );

        var_num_.push_back( aux.add_variable("f_x", order, libMesh::LAGRANGE) );
        var_num_.push_back( aux.add_variable("f_y", order, libMesh::LAGRANGE) );

        if(dim > 2)
            var_num_.push_back( aux.add_variable("f_z", order, libMesh::LAGRANGE) );

        var_num_.push_back( aux.add_variable("fext_x", order, libMesh::LAGRANGE) );
        var_num_.push_back( aux.add_variable("fext_y", order, libMesh::LAGRANGE) );

        if(dim > 2)
            var_num_.push_back( aux.add_variable("fext_z", order, libMesh::LAGRANGE) );


        var_num_.push_back( aux.add_variable("n_x", order, libMesh::LAGRANGE) );
        var_num_.push_back( aux.add_variable("n_y", order, libMesh::LAGRANGE) );

        if(dim > 2)
            var_num_.push_back( aux.add_variable("n_z", order, libMesh::LAGRANGE) );

        var_num_.push_back( aux.add_variable("wear_disp_x", order, libMesh::LAGRANGE) );
        var_num_.push_back( aux.add_variable("wear_disp_y", order, libMesh::LAGRANGE) );

        if(dim > 2)
            var_num_.push_back( aux.add_variable("wear_disp_z", order, libMesh::LAGRANGE) );


        var_num_.push_back( aux.add_variable("normalstress", order, libMesh::LAGRANGE) );
        var_num_.push_back( aux.add_variable("slidingdistance", order, libMesh::LAGRANGE) );
        var_num_.push_back( aux.add_variable("wear", order, libMesh::LAGRANGE) );
        var_num_.push_back( aux.add_variable("is_contact", order, libMesh::LAGRANGE) );
        var_num_.push_back( aux.add_variable("gap", order, libMesh::LAGRANGE) );

        aux.init();
        aux.update();

        wear_ = local_zeros(dof_map_main.n_local_dofs());
    }

    void ContactSystem::update_wear(
        const double dt,
        const double wear_coefficient,
        const UVector &sliding_distance,
        const UVector &normal_stress,
        UVector &wear
        )
    {
        wear += (dt * wear_coefficient) * abs(e_mul(sliding_distance, normal_stress));
    }

    void ContactSystem::update(
        const MechanicsState &state,
        const ContactT &contact,
        const double dt)
    {
        using libMesh::MeshBase;
        auto &aux = equation_systems_->get_system("aux");
        const auto &mesh = equation_systems_->get_mesh();
        const int dim = mesh.mesh_dimension();

        UVector normal_stress = local_zeros(local_size(state.displacement));
        UVector sliding_distance = local_zeros(local_size(state.displacement));

            // UVector stress = local_zeros(local_size(state.displacement));

        if(contact.initialized()) {
            // normal_stress = contact.orthogonal_trafo * state.stress;

            contact.apply_orthogonal_trafo(state.stress, normal_stress);

            // UVector tangential_velocity = contact.orthogonal_trafo * state.velocity;

            UVector tangential_velocity;
            contact.apply_orthogonal_trafo(state.velocity, tangential_velocity);
            sliding_distance = local_zeros(local_size(state.displacement));

            {
                Read<UVector>   r_v(tangential_velocity);
                Write<UVector> w_s(sliding_distance);

                Range r = range(tangential_velocity);

                for(auto i = r.begin(); i < r.end(); i += dim) {

                    double dist = 0.;
                    for(uint d = 1; d < dim; ++d) {
                        auto val = tangential_velocity.get(i + d);
                        dist += val*val;
                    }

                    dist = std::sqrt(dist);
                    sliding_distance.set(i, dist);
                }
            }

            sliding_distance = e_mul(contact.is_contact_node(), sliding_distance);
                // normal_stress = e_mul(contact.is_contact_node, contact.orthogonal_trafo * stress);

                //override normal stress
            {
                normal_stress *= 0.;

                Read<UVector> r_s(state.stress), r_n(contact.normals());
                Write<UVector> w_n(normal_stress);

                Range r = range(state.stress);

                for(auto i = r.begin(); i != r.end(); i+= dim) {
                    double ns = 0.;
                    for(unsigned int d = 0; d < dim; ++d) {
                        ns += state.stress.get(i+d) * contact.normals().get(i+d);
                    }

                    normal_stress.set(i, ns);
                }
            }

            update_wear(dt, wear_coefficient_, sliding_distance, normal_stress, wear_);
            total_wear_.push_back(sum(wear_));
        }

        {
            Read<UVector> r_d(state.displacement), r_v(state.velocity);
            Read<UVector> r_f(state.internal_force), r_ef(state.external_force), r_ns(normal_stress);
            Read<UVector> r_c(contact.is_contact_node());

            auto nd 	= mesh.local_nodes_begin();
            auto nd_end = mesh.local_nodes_end();

            for (; nd != nd_end; ++nd) {
                const libMesh::Node * node = *nd;

                for (unsigned int d = 0; d < dim; ++d) {
                    unsigned int source_dof = node->dof_number(main_system_num_, d, 0);

                    auto dest_dof_disp 	= node->dof_number(aux.number(), var_num_[d], 0);
                    auto dest_dof_vel 	= node->dof_number(aux.number(), var_num_[dim + d], 0);
                    auto dest_dof_force = node->dof_number(aux.number(), var_num_[2 * dim + d], 0);
                    auto dest_dof_ext_force = node->dof_number(aux.number(), var_num_[3 * dim + d], 0);

                    aux.solution->set(dest_dof_disp,  state.displacement_increment.get(source_dof));
                    aux.solution->set(dest_dof_vel,   state.velocity.get(source_dof));
                        // aux.solution->set(dest_dof_force, state.internal_force.get(source_dof));

                    if(contact.initialized()) {
                        aux.solution->set(dest_dof_force, state.stress.get(source_dof));
                    }

                    aux.solution->set(dest_dof_ext_force, state.external_force.get(source_dof));

                    auto dest_dof_n = node->dof_number(aux.number(), var_num_[4 * dim + d], 0);

                    if(contact.initialized()) {
                        aux.solution->set(dest_dof_n, contact.normals().get(source_dof));
                    } else {
                        aux.solution->set(dest_dof_n, 0.);
                    }

                    // auto dest_dof_wear_disp = node->dof_number(aux.number(), var_num_[5 * dim + d], 0);
                    // aux.solution->set(dest_dof_wear_disp, wear_induced_displacement.get(source_dof));
                }

                unsigned int source_dof = node->dof_number(main_system_num_, 0, 0);
                unsigned int dest_dof_normal_stress = node->dof_number(aux.number(), var_num_[6 * dim], 0);
                aux.solution->set(dest_dof_normal_stress, normal_stress.get(source_dof));

                unsigned int dest_dof_sliding_dist = node->dof_number(aux.number(), var_num_[6 * dim + 1], 0);
                aux.solution->set(dest_dof_sliding_dist, sliding_distance.get(source_dof));

                unsigned int dest_dof_wear = node->dof_number(aux.number(), var_num_[6 * dim + 2], 0);
                aux.solution->set(dest_dof_wear, wear_.get(source_dof));

                unsigned int dest_dof_is_contact = node->dof_number(aux.number(), var_num_[6 * dim + 3], 0);
                unsigned int dest_dof_is_gap = node->dof_number(aux.number(), var_num_[6 * dim + 4], 0);

                if(contact.initialized()) {
                    aux.solution->set(dest_dof_is_contact, contact.is_contact_node().get(source_dof));
                    aux.solution->set(dest_dof_is_gap, contact.gap().get(source_dof));
                } else {
                    aux.solution->set(dest_dof_is_contact, 0.);
                    aux.solution->set(dest_dof_is_gap, 0.);
                }
            }
        }

        aux.solution->close();
    }
}

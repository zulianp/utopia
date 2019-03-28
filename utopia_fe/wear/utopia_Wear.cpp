#include "utopia_Wear.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_NormalTangentialCoordinateSystem.hpp"

#include "moonolith_communicator.hpp"
#include "moonolith_describe.hpp"
#include "moonolith_synched_describable.hpp"
#include "utopia_Contact.hpp"


namespace utopia {

    void apply_displacement(
        const UVector &displacement_increment,
        const libMesh::DofMap &dof_map,
        libMesh::MeshBase &mesh)
    {
        //FIXME
        int sys_num = 0;

        moonolith::Communicator comm(mesh.comm().get());

        if(!comm.is_alone()) {
            auto r = range(displacement_increment);
            Read<UVector> r_d(displacement_increment);

            auto m_begin = mesh.active_local_elements_begin();
            auto m_end   = mesh.active_local_elements_end();

            std::vector<PetscInt> idx;
            std::set<PetscInt> unique_idx;
            std::map<libMesh::dof_id_type, double> idx_to_value;
            std::vector<libMesh::dof_id_type> dof_indices;

            for(auto m_it = m_begin; m_it != m_end; ++m_it) {
                dof_map.dof_indices(*m_it, dof_indices);
                for(auto dof_id : dof_indices) {
                    if(r.inside(dof_id)) {
                        idx_to_value[dof_id] = displacement_increment.get(dof_id);
                    } else {
                        unique_idx.insert(dof_id);
                    }
                }
            }

            idx.insert(idx.end(), unique_idx.begin(), unique_idx.end());
            UVector out = displacement_increment.select(idx);
            {
                Read<UVector> r_out(out);
                auto range_out = range(out);

                for(std::size_t i = 0; i < idx.size(); ++i) {
                    idx_to_value[idx[i]] = out.get(range_out.begin() + i);
                }
            }

            for(auto m_it = m_begin; m_it != m_end; ++m_it) {
                auto &e = **m_it;
                for(int i = 0; i < e.n_nodes(); ++i) {
                    auto &node = e.node_ref(i);

                    for(unsigned int c = 0; c < mesh.mesh_dimension(); ++c) {
                        const int dof_id = node.dof_number(sys_num, c, 0);
                        assert(idx_to_value.find(dof_id) != idx_to_value.end());
                        double &val = idx_to_value[dof_id];
                        node(c) += val;
                        val = 0.;
                    }
                }
            }

        } else {

            Read<UVector> r_d(displacement_increment);

            auto m_it  = mesh.local_nodes_begin();
            auto m_end = mesh.local_nodes_end();


            for(; m_it != m_end; ++m_it) {
                for(unsigned int c = 0; c < mesh.mesh_dimension(); ++c) {
                    const int dof_id = (*m_it)->dof_number(sys_num, c, 0);
                    (**m_it)(c) += displacement_increment.get(dof_id);
                }
            }
        }
    }

    Wear::Wear()
    : wear_coefficient(7e-3), extrapolation_factor(10.)
    {

    }

    void Wear::compute_displacement(
        ProductFunctionSpace<LibMeshFunctionSpace> &V,
        const std::vector<int> &boundary_tags,
        UVector &wear_induced_displacement
        )
    {
        libMesh::MeshBase &mesh  = V[0].mesh();
        libMesh::DofMap &dof_map = V[0].dof_map();
        auto dim = mesh.mesh_dimension();

        wear_induced_displacement = local_zeros(dof_map.n_local_dofs());

        //compute surface normals
        assemble_normal_tangential_transformation(mesh, dof_map, boundary_tags, is_normal_component, normals, trafo);

        //compute surface displacement
        // wear_induced_displacement = local_zeros(local_size(wear));
        {
            auto r = range(wear);
            Write<UVector> w_w(wear_induced_displacement);
            Read<UVector> r_w(wear), r_n(normals), r_i(is_normal_component);

            for(auto i = r.begin(); i < r.end(); i += dim) {
                if(is_normal_component.get(i) > 0) {
                    for(unsigned int d = 0; d < dim; ++d) {
                        assert(wear.get(i) >= -1e-16);
                        wear_induced_displacement.set(i + d, normals.get(i + d) * (-extrapolation_factor) * wear.get(i));
                    }
                }
            }
        }

    }

    void Wear::mesh_displacement(
        ProductFunctionSpace<LibMeshFunctionSpace> &V,
        const std::vector<int> &boundary_tags,
        UVector &warped_displacement)
    {
        libMesh::MeshBase &mesh = V[0].mesh();
        libMesh::DofMap &dof_map = V[0].dof_map();
        auto dim = mesh.mesh_dimension();

        compute_displacement(V, boundary_tags, wear_induced_displacement);

        // the interior using linear elasticity or laplacian
        // use dirichlet conditions on the whole boundary and warp
        auto u = trial(V);
        auto v = test(V);

        USparseMatrix lapl_mat;
        auto lapl = inner(grad(u), grad(v)) * dX;
        assemble(lapl, lapl_mat);
        warped_displacement = local_zeros(local_size(wear_induced_displacement));

        //FIXME warped_displacement passed as dummy
        set_identity_at_constraint_rows(dof_map, lapl_mat);

        // KSPSolver<USparseMatrix, UVector> solver;
        Factorization<USparseMatrix, UVector> solver;
        if(!solver.solve(lapl_mat, wear_induced_displacement, warped_displacement)) {
            std::cerr << "[Warning] harmonic map did not work" << std::endl;
            warped_displacement = wear_induced_displacement;
        }

        //displace mesh
        // apply_displacement(warped_displacement, dof_map, mesh);
        convert(warped_displacement, *V[0].equation_system().solution);
        // convert(wear_induced_displacement, *V[0].equation_system().solution);
        V[0].equation_system().solution->close();

        // double wear_magnitude = norm2(wear_induced_displacement);
        // std::cout << "wear_magnitude: " << wear_magnitude << std::endl;

        // double param_magn = norm2(warped_displacement);
        // std::cout << "param_magn: " << param_magn << std::endl;
    }

    void Wear::init_aux_system(
        libMesh::EquationSystems &es,
        libMesh::Order order)
    {
        //init aux system for plotting
        auto &aux = es.add_system<libMesh::LinearImplicitSystem>("aux");

        var_num_aux.push_back( aux.add_variable("inc_x", libMesh::Order(order), libMesh::LAGRANGE) );
        var_num_aux.push_back( aux.add_variable("inc_y", libMesh::Order(order), libMesh::LAGRANGE) );

        const int dim = es.get_mesh().mesh_dimension();

        if(dim > 2)
            var_num_aux.push_back( aux.add_variable("inc_z", libMesh::Order(order), libMesh::LAGRANGE) );

        var_num_aux.push_back( aux.add_variable("vel_x", libMesh::Order(order), libMesh::LAGRANGE) );
        var_num_aux.push_back( aux.add_variable("vel_y", libMesh::Order(order), libMesh::LAGRANGE) );

        if(dim > 2)
            var_num_aux.push_back( aux.add_variable("vel_z", libMesh::Order(order), libMesh::LAGRANGE) );

        var_num_aux.push_back( aux.add_variable("f_x", libMesh::Order(order), libMesh::LAGRANGE) );
        var_num_aux.push_back( aux.add_variable("f_y", libMesh::Order(order), libMesh::LAGRANGE) );

        if(dim > 2)
            var_num_aux.push_back( aux.add_variable("f_z", libMesh::Order(order), libMesh::LAGRANGE) );

        var_num_aux.push_back( aux.add_variable("fext_x", libMesh::Order(order), libMesh::LAGRANGE) );
        var_num_aux.push_back( aux.add_variable("fext_y", libMesh::Order(order), libMesh::LAGRANGE) );

        if(dim > 2)
            var_num_aux.push_back( aux.add_variable("fext_z", libMesh::Order(order), libMesh::LAGRANGE) );


        var_num_aux.push_back( aux.add_variable("n_x", libMesh::Order(order), libMesh::LAGRANGE) );
        var_num_aux.push_back( aux.add_variable("n_y", libMesh::Order(order), libMesh::LAGRANGE) );

        if(dim > 2)
            var_num_aux.push_back( aux.add_variable("n_z", libMesh::Order(order), libMesh::LAGRANGE) );

        var_num_aux.push_back( aux.add_variable("wear_disp_x", libMesh::Order(order), libMesh::LAGRANGE) );
        var_num_aux.push_back( aux.add_variable("wear_disp_y", libMesh::Order(order), libMesh::LAGRANGE) );

        if(dim > 2)
            var_num_aux.push_back( aux.add_variable("wear_disp_z", libMesh::Order(order), libMesh::LAGRANGE) );


        var_num_aux.push_back( aux.add_variable("normalstress", libMesh::Order(order), libMesh::LAGRANGE) );
        var_num_aux.push_back( aux.add_variable("slidingdistance", libMesh::Order(order), libMesh::LAGRANGE) );
        var_num_aux.push_back( aux.add_variable("wear", libMesh::Order(order), libMesh::LAGRANGE) );
        var_num_aux.push_back( aux.add_variable("is_contact", libMesh::Order(order), libMesh::LAGRANGE) );
        var_num_aux.push_back( aux.add_variable("gap", libMesh::Order(order), libMesh::LAGRANGE) );

        aux.init();
        aux.update();
    }

    void Wear::update_aux_system(
        const int main_system_number,
        const MechanicsState &state,
        const Contact &contact,
        const double dt,
        libMesh::EquationSystems &es)
    {
        using libMesh::MeshBase;
        auto &main = es.get_system(main_system_number);
        auto &aux = es.get_system("aux");
        const auto &mesh = es.get_mesh();
        const int dim = mesh.mesh_dimension();

        UVector normal_stress = local_zeros(local_size(state.displacement));
        UVector sliding_distance = local_zeros(local_size(state.displacement));

        if(empty(wear_induced_displacement)) {
            wear_induced_displacement = local_zeros(local_size(state.displacement));
        }

        // UVector stress = local_zeros(local_size(state.displacement));

        if(contact.initialized) {
            // stress = e_mul(mech_ctx.inverse_mass_vector, (state.external_force - state.internal_force));
            // stress = (state.external_force - state.internal_force);
            normal_stress = contact.orthogonal_trafo * state.stress;

            UVector tangential_velocity = contact.orthogonal_trafo * state.velocity;
            sliding_distance = local_zeros(local_size(state.velocity));

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

            sliding_distance = e_mul(contact.is_contact_node, sliding_distance);
            // normal_stress = e_mul(contact.is_contact_node, contact.orthogonal_trafo * stress);

            //override normal stress
            {
                normal_stress *= 0.;

                Read<UVector> r_s(state.stress), r_n(contact.normals);
                Write<UVector> w_n(normal_stress);

                Range r = range(state.stress);

                for(auto i = r.begin(); i != r.end(); i+= dim) {
                    double ns = 0.;
                    for(unsigned int d = 0; d < dim; ++d) {
                        ns += state.stress.get(i+d) * contact.normals.get(i+d);
                    }

                    normal_stress.set(i, ns);
                }
            }

            const double mag_normal_stress = norm2(normal_stress);
            const double mag_sliding_dist  = norm2(sliding_distance);

            std::cout << "mag_normal_stress: " << mag_normal_stress << std::endl;
            std::cout << "mag_sliding_dist:  " << mag_sliding_dist  << std::endl;

            update(dt, sliding_distance, normal_stress);
            total_wear.push_back(sum(wear));
        }

        {
            Read<UVector> r_d(state.displacement), r_v(state.velocity);
            Read<UVector> r_f(state.internal_force), r_ef(state.external_force), r_ns(normal_stress);
            Read<UVector> r_c(contact.is_contact_node);

            auto nd 	= mesh.local_nodes_begin();
            auto nd_end = mesh.local_nodes_end();

            for (; nd != nd_end; ++nd) {
                const libMesh::Node * node = *nd;

                for (unsigned int d = 0; d < dim; ++d) {
                    unsigned int source_dof = node->dof_number(main_system_number, d, 0);

                    auto dest_dof_disp 	= node->dof_number(aux.number(), var_num_aux[d], 0);
                    auto dest_dof_vel 	= node->dof_number(aux.number(), var_num_aux[dim + d], 0);
                    auto dest_dof_force = node->dof_number(aux.number(), var_num_aux[2 * dim + d], 0);
                    auto dest_dof_ext_force = node->dof_number(aux.number(), var_num_aux[3 * dim + d], 0);

                    aux.solution->set(dest_dof_disp,  state.displacement_increment.get(source_dof));
                    aux.solution->set(dest_dof_vel,   state.velocity.get(source_dof));
                    // aux.solution->set(dest_dof_force, state.internal_force.get(source_dof));

                    if(contact.initialized) {
                        aux.solution->set(dest_dof_force, state.stress.get(source_dof));
                    }

                    aux.solution->set(dest_dof_ext_force, state.external_force.get(source_dof));

                    auto dest_dof_n = node->dof_number(aux.number(), var_num_aux[4 * dim + d], 0);

                    if(contact.initialized) {
                        aux.solution->set(dest_dof_n, contact.normals.get(source_dof));
                    } else {
                        aux.solution->set(dest_dof_n, 0.);
                    }

                    auto dest_dof_wear_disp = node->dof_number(aux.number(), var_num_aux[5 * dim + d], 0);
                    aux.solution->set(dest_dof_wear_disp, wear_induced_displacement.get(source_dof));
                }

                unsigned int source_dof = node->dof_number(main_system_number, 0, 0);
                unsigned int dest_dof_normal_stress = node->dof_number(aux.number(), var_num_aux[6 * dim], 0);
                aux.solution->set(dest_dof_normal_stress, normal_stress.get(source_dof));

                unsigned int dest_dof_sliding_dist = node->dof_number(aux.number(), var_num_aux[6 * dim + 1], 0);
                aux.solution->set(dest_dof_sliding_dist, sliding_distance.get(source_dof));

                unsigned int dest_dof_wear = node->dof_number(aux.number(), var_num_aux[6 * dim + 2], 0);
                aux.solution->set(dest_dof_wear, wear.get(source_dof));

                unsigned int dest_dof_is_contact = node->dof_number(aux.number(), var_num_aux[6 * dim + 3], 0);
                unsigned int dest_dof_is_gap = node->dof_number(aux.number(), var_num_aux[6 * dim + 4], 0);

                if(contact.initialized) {
                    aux.solution->set(dest_dof_is_contact, contact.is_contact_node.get(source_dof));
                    aux.solution->set(dest_dof_is_gap, contact.gap.get(source_dof));
                } else {
                    aux.solution->set(dest_dof_is_contact, 0.);
                    aux.solution->set(dest_dof_is_gap, 0.);
                }
            }
        }

        aux.solution->close();
    }
}

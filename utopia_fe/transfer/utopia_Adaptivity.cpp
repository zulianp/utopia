
#include "utopia_Adaptivity.hpp"
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_elem.h"
#include "libmesh/parallel_node.h"
#include "libmesh/remote_elem.h"
// #include "libmesh/parallel_sync.h"
#include "libmesh/coupling_matrix.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/numeric_vector.h"

#include "libmesh/boundary_info.h"

#include "utopia_LibMeshBackend.hpp"

namespace utopia {

    void Adaptivity::compute_all_constraints(const libMesh::MeshBase &mesh,
                                             const libMesh::DofMap &dof_map,
                                             libMesh::DofConstraints &constraints) {
        using uint = unsigned int;

        std::vector<int> index;

        libMesh::MeshBase &mesh_copy = const_cast<libMesh::MeshBase &>(mesh);

        libMesh::DofMap &dof_copy = const_cast<libMesh::DofMap &>(dof_map);

        uint vars = dof_map.n_variables();

        uint var_num = 0;

        // for(uint var_num = 0; var_num < vars; ++var_num) {

        libMesh::FEType fe_type = dof_map.variable_type(var_num);

        fe_type.order = static_cast<libMesh::Order>(fe_type.order);

        if (fe_type.order > 0) {
            compute_boundary_nodes(mesh, dof_copy, 0, 0, index);
        }
        // }

        auto el = mesh.active_local_elements_begin();

        const auto end_el = mesh.active_local_elements_end();

        for (; el != end_el; ++el) {
            const auto *elem = *el;

            auto *ele = *el;

            libmesh_assert(mesh.is_prepared());

            uint n_vars = dof_map.n_variables();

            uint var_num = 0;

            // for(uint var_num = 0; var_num < n_vars; ++var_num) {

            libMesh::FEType fe_type = dof_map.variable_type(var_num);

            fe_type.order = static_cast<libMesh::Order>(fe_type.order);

            if (fe_type.order > 0) {
                compute_constraints(constraints, dof_map, var_num, elem, mesh.mesh_dimension());
            }
            // }
        }

        if (fe_type.order > 0) {
            if (mesh.mesh_dimension() == 3) dof_copy.process_constraints(mesh_copy);

            process_constraints(mesh_copy, dof_copy, constraints, index);
        }

        std::cout << "--------------------------------------------------\n";
        std::cout << "[Adaptivity::compute_all_constraints] n_constraints: " << constraints.size() << std::endl;
        std::cout << "--------------------------------------------------\n";
    }

    void Adaptivity::assemble_constraint(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map) {
        auto el = mesh.active_local_elements_begin();

        const auto end_el = mesh.active_local_elements_end();

        dof_constraints_.clear();

        std::vector<int> index;

        libMesh::MeshBase &mesh_copy = const_cast<libMesh::MeshBase &>(mesh);

        libMesh::DofMap &dof_copy = const_cast<libMesh::DofMap &>(dof_map);

        uint var_num = 0;

        // for(uint var_num = 0; var_num < vars; ++var_num) {

        libMesh::FEType fe_type = dof_map.variable_type(var_num);

        fe_type.order = static_cast<libMesh::Order>(fe_type.order);

        if (fe_type.order > 0) {
            compute_boundary_nodes(mesh, dof_copy, 0, 0, index);
        }

        if (fe_type.order > 0) {
            for (; el != end_el; ++el) {
                const auto *elem = *el;

                auto *ele = *el;

                libmesh_assert(mesh.is_prepared());

                uint var_num = 0;

                // for (unsigned int var_num=0; var_num<dof_map.n_variables();
                //      ++var_num) {
                compute_constraints(dof_constraints_, dof_map, var_num, elem, mesh.mesh_dimension());
                //ß}
            }

            if (mesh.mesh_dimension() == 3) dof_copy.process_constraints(mesh_copy);

            Adaptivity::process_constraints(mesh_copy, dof_copy, dof_constraints_, index);
        }

        std::cout << "--------------------------------------------------\n";
        std::cout << "[Adaptivity::assemble_constraint] n_constraints: " << dof_constraints_.size() << std::endl;
        std::cout << "--------------------------------------------------\n";
    }

    void Adaptivity::constraint_matrix(const LibMeshFunctionSpace &V, USparseMatrix &M, USparseMatrix &S) {
        const auto &mesh = V.mesh();
        unsigned int var_num = V.subspace_id();

        const auto &dof_map = V.dof_map();
        constraint_matrix(mesh, dof_map, M, S);
    }

    void Adaptivity::constraint_matrix(const libMesh::MeshBase &mesh,
                                       const libMesh::DofMap &dof_map,
                                       USparseMatrix &M,
                                       USparseMatrix &S) {
        std::vector<SizeType> index;

        assemble_constraint(mesh, dof_map);

        bool called_recursively = false;

        M = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), 30);

        S = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), 30);

        std::vector<libMesh::dof_id_type> elem_dofs;

        std::vector<libMesh::dof_id_type> I(1, 0), J(1, 0);

        std::vector<double> V(1, 0);

        auto el = mesh.active_local_elements_begin();

        const auto end_el = mesh.active_local_elements_end();

        libMesh::DenseMatrix<libMesh::Number> C;

        Write<USparseMatrix> w(M, utopia::GLOBAL_INSERT), w_2(S, utopia::GLOBAL_INSERT);

        for (; el != end_el; ++el) {
            const auto *elem = *el;

            auto *ele = *el;

            dof_map.dof_indices(elem, elem_dofs);

            std::set<libMesh::dof_id_type> dof_set;

            bool we_have_constraints = false;

            for (const auto &dof : elem_dofs) {
                if (dof_map.is_constrained_dof(dof)) {
                    we_have_constraints = true;

                    auto pos = dof_constraints_.find(dof);

                    if (pos == dof_constraints_.end()) {
                        M.c_set(dof, dof, 1.0);

                        continue;
                    }

                    const auto &constraint_row = pos->second;

                    for (const auto &item : constraint_row) dof_set.insert(item.first);
                }
            }

            for (const auto &dof : elem_dofs) dof_set.erase(dof);

            if (!dof_set.empty() || !called_recursively) {
                const unsigned int old_size = static_cast<unsigned int>(elem_dofs.size());

                elem_dofs.insert(elem_dofs.end(), dof_set.begin(), dof_set.end());

                C.resize(old_size, static_cast<unsigned int>(elem_dofs.size()));

                for (unsigned int i = 0; i != old_size; i++)
                    if (dof_map.is_constrained_dof(elem_dofs[i])) {
                        // I[0] = elem_dofs[i];
                        // J[0] = elem_dofs[i];
                        // V[0] = 0.0;

                        // //S.set (elem_dofs[i],elem_dofs[i],0.0);
                        // S.set_matrix(I, J, V);
                        S.c_set(elem_dofs[i], elem_dofs[i], 0.0);

                        auto pos = dof_constraints_.find(elem_dofs[i]);

                        if (pos == dof_constraints_.end()) continue;

                        const auto &constraint_row = pos->second;

                        for (const auto &item : constraint_row) {
                            for (unsigned int j = 0, n_elem_dofs = static_cast<unsigned int>(elem_dofs.size());
                                 j != n_elem_dofs;
                                 j++) {
                                if (elem_dofs[j] == item.first) {
                                    C(i, j) = item.second;
                                    M.c_set(elem_dofs[i], elem_dofs[j], item.second);
                                    // M.set (elem_dofs[i],elem_dofs[j],item.second);

                                    S.c_set(elem_dofs[i], elem_dofs[j], item.second);
                                    // S.set (elem_dofs[i],elem_dofs[j],item.second);
                                }
                            }
                        }
                    }

                    else {
                        C(i, i) = 1.;

                        M.c_set(elem_dofs[i], elem_dofs[i], 1.0);
                        // M.set (elem_dofs[i],elem_dofs[i],1.0);
                    }
            }
        }
    }

    void Adaptivity::compute_constraints(libMesh::DofConstraints &constraints,
                                         const libMesh::DofMap &dof_map,
                                         const unsigned int var_num,
                                         const libMesh::Elem *elem,
                                         const unsigned mesh_dim) {

     #ifdef LIBMESH_VERSION_LESS_THAN(1, 4, 0)
     #else
        // Only constrain elements in 2,3D.
        if (mesh_dim == 1) return;

        // std::cout<<"lagrange_compute_constraints libmesh mio
        // prima:"<<constraints.size()<<std::endl;
        libmesh_assert(elem);

        // Only constrain active and ancestor elements
        if (elem->subactive())  // if the element is subactive (i.e. has no active
                                // descendants)
            return;

        libMesh::FEType fe_type = dof_map.variable_type(var_num);
        fe_type.order = static_cast<libMesh::Order>(fe_type.order + elem->p_level());

        std::vector<libMesh::dof_id_type> my_dof_indices, parent_dof_indices;
        libMesh::UniquePtr<const libMesh::Elem> my_side, parent_side;

        for (auto s : elem->side_index_range()) {
            if (elem->neighbor_ptr(s) != nullptr && elem->neighbor_ptr(s) != libMesh::remote_elem) {
                if (elem->neighbor_ptr(s)->level() < elem->level())

                {
                    const auto *parent = elem->parent();
                    libmesh_assert(parent);

                    elem->build_side_ptr(my_side, s);

                    parent->build_side_ptr(parent_side, s);

                    my_dof_indices.reserve(my_side->n_nodes());

                    parent_dof_indices.reserve(parent_side->n_nodes());

                    dof_map.dof_indices(my_side.get(), my_dof_indices, var_num);

                    dof_map.dof_indices(parent_side.get(), parent_dof_indices, var_num);

                    const unsigned int n_side_dofs =
                        libMesh::FEInterface::n_dofs(mesh_dim - 1, fe_type, my_side->type());

                    const unsigned int n_parent_side_dofs =
                        libMesh::FEInterface::n_dofs(mesh_dim - 1, fe_type, parent_side->type());

                    for (unsigned int my_dof = 0; my_dof != n_side_dofs; my_dof++) {
                        libmesh_assert_less(my_dof, my_side->n_nodes());
                        assert(my_dof < n_side_dofs);

                        const libMesh::dof_id_type my_dof_g = my_dof_indices[my_dof];

                        bool self_constraint = false;

                        for (unsigned int their_dof = 0; their_dof != n_parent_side_dofs; their_dof++) {
                            libmesh_assert_less(their_dof, parent_side->n_nodes());

                            const libMesh::dof_id_type their_dof_g = parent_dof_indices[their_dof];

                            if (their_dof_g == my_dof_g) {
                                self_constraint = true;
                                break;
                            }
                        }

                        if (self_constraint) continue;

                        libMesh::DofConstraintRow *constraint_row;

                        constraint_row = &(constraints[my_dof_g]);

                        const libMesh::Point &support_point = my_side->point(my_dof);

                        const libMesh::Point mapped_point =
                            libMesh::FEInterface::inverse_map(mesh_dim - 1, fe_type, parent_side.get(), support_point);

                        for (unsigned int their_dof = 0; their_dof != n_parent_side_dofs; their_dof++) {
                            libmesh_assert_less(their_dof, parent_side->n_nodes());

                            const libMesh::dof_id_type their_dof_g = parent_dof_indices[their_dof];

                            const double their_dof_value = libMesh::FEInterface::shape(
                                mesh_dim - 1, fe_type, parent_side->type(), their_dof, mapped_point);
                            if ((std::abs(their_dof_value) > 1.e-5) && (std::abs(their_dof_value) < .999)) {
                                // auto check = std::find(index.begin(), index.end(), their_dof_g)
                                // != index.end();

                                // if(!check)
                                constraint_row->insert(std::make_pair(their_dof_g, their_dof_value));
                            }
                        }
                    }
                }
            }
        }
     #endif
    }

    void Adaptivity::process_constraints(libMesh::MeshBase &mesh,
                                         libMesh::DofMap &dof_map,
                                         libMesh::DofConstraints &_dof_constraints,
                                         std::vector<int> &index) {
        using namespace libMesh;

        // dof_map.allgather_recursive_constraints(mesh);

        // auto _primal_constraint_values = dof_map.get_primal_constraint_values();

        // std::vector<int> index;

        typedef std::set<dof_id_type> RCSet;

        RCSet unexpanded_set;

        for (const auto &i : _dof_constraints) unexpanded_set.insert(i.first);

        while (!unexpanded_set.empty())
            for (RCSet::iterator i = unexpanded_set.begin(); i != unexpanded_set.end();
                 /* nothing */) {
                // If the DOF is constrained
                DofConstraints::iterator pos = _dof_constraints.find(*i);

                libmesh_assert(pos != _dof_constraints.end());

                DofConstraintRow &constraint_row = pos->second;

                std::vector<dof_id_type> constraints_to_expand;

                for (const auto &item : constraint_row) {
                    if (item.first != *i && dof_map.is_constrained_dof(item.first)) {
                        bool check = true;

                        for (auto it = index.begin(); it < index.end(); ++it) {
                            int b_id = *it;

                            if (b_id == item.first) {
                                check = false;
                            }
                        }

                        if (check == true) {
                            unexpanded_set.insert(item.first);

                            constraints_to_expand.push_back(item.first);
                        }
                    }
                }

                for (const auto &expandable : constraints_to_expand) {
                    const Real this_coef = constraint_row[expandable];

                    DofConstraints::const_iterator subpos = _dof_constraints.find(expandable);

                    libmesh_assert(subpos != _dof_constraints.end());

                    if (subpos == _dof_constraints.end()) return;

                    const DofConstraintRow &subconstraint_row = subpos->second;

                    for (const auto &item : subconstraint_row) {
                        // Assert that the constraint does not form a cycle.
                        libmesh_assert(item.first != expandable);
                        constraint_row[item.first] += item.second * this_coef;
                    }

                    constraint_row.erase(expandable);
                }

                if (constraints_to_expand.empty())
                    i = unexpanded_set.erase(i);
                else
                    ++i;
            }
    }

    void Adaptivity::compute_boundary_nodes(const libMesh::MeshBase &mesh,
                                            libMesh::DofMap &dof_map,
                                            unsigned int sys_number,
                                            unsigned int var_number,
                                            std::vector<int> &index) {
        // std::cout << "Adaptivity::compute_boundary_nodes::Begin " << std::endl;
      #ifdef LIBMESH_VERSION_LESS_THAN(1, 4, 0)
      #else
        auto on_boundary = libMesh::MeshTools::find_boundary_nodes(mesh);

        std::vector<int> dirichlet_id, index_local;

        index_local.clear();

        dirichlet_id.clear();

        index.clear();

        if (mesh.mesh_dimension() < 3) {
            libMesh::MeshBase::const_node_iterator it = mesh.local_nodes_begin();
            ;

            const libMesh::MeshBase::const_node_iterator end_it = mesh.local_nodes_end();

            for (; it != end_it; ++it) {
                // const libMesh::Elem * ele = *it;

                const libMesh::Node *node = *it;

                // for(int kk=0; kk<ele->n_sides(); kk++) {

                // auto neigh = ele->neighbor_ptr(kk);

                // if (neigh != libMesh::remote_elem &&
                // utopia::boundary_id(mesh.get_boundary_info(), ele, kk)>0)
                //{
                // auto side = ele->build_side_ptr(kk);

                index_local.clear();

                // for (int ll=0; ll<ele->n_nodes(); ll++)
                // {

                // const libMesh::Node * node = ele->node_ptr(ll);

                const libMesh::dof_id_type node_dof = node->dof_number(sys_number, var_number, 0);

                if (on_boundary.count(node->id()) && dof_map.is_constrained_dof(node_dof)) {
                    index_local.push_back(node_dof);

                    index.push_back(node_dof);

                    // std::cout<<node_dof<<std::endl;
                }

                // }

                // if(index_local.size()==side->n_nodes())
                // {
                // index.insert(index.end(), index_local.begin(), index_local.end());
                //}
                //}
                //}
            }
        }

        else {
            {
                libMesh::DofConstraintValueMap &rhs_values = dof_map.get_primal_constraint_values();

                libMesh::MeshBase::const_element_iterator it = mesh.active_elements_begin();

                const libMesh::MeshBase::const_element_iterator end_it = mesh.active_elements_end();

                std::vector<libMesh::dof_id_type> my_dof_indices, parent_dof_indices;

                libMesh::UniquePtr<const libMesh::Elem> my_side, parent_side_0;

                libMesh::FEType fe_type = dof_map.variable_type(0);

                for (; it != end_it; ++it) {
                    const libMesh::Elem *ele = *it;

                    const auto *ele_parent_0 = ele->top_parent();

                    for (int jj = 0; jj < ele_parent_0->n_sides(); jj++) {
                        libmesh_assert(ele_parent_0);

                        auto parent_side_0 = ele_parent_0->build_side_ptr(jj);

                        index_local.clear();

                        for (int ll = 0; ll < parent_side_0->n_nodes(); ll++) {
                            const libMesh::Node *node_0 = parent_side_0->node_ptr(ll);

                            const libMesh::dof_id_type node_dof_0 = node_0->dof_number(sys_number, var_number, 0);

                            if (dof_map.is_constrained_dof(node_dof_0)) {
                                index_local.push_back(node_dof_0);

                                auto valpos = rhs_values.find(node_dof_0);

                                index.push_back(node_dof_0);
                            }
                        }

                        if (index_local.size() == parent_side_0->n_nodes()) {
                            auto bc_id = utopia::boundary_id(mesh.get_boundary_info(), ele_parent_0, jj);

                            auto check =
                                (std::find(dirichlet_id.begin(), dirichlet_id.end(), bc_id) != dirichlet_id.end());

                            if (!check) dirichlet_id.push_back(bc_id);
                        }
                    }
                }
            }

            libMesh::MeshBase::const_element_iterator it = mesh.active_elements_begin();

            const libMesh::MeshBase::const_element_iterator end_it = mesh.active_elements_end();

            for (; it != end_it; ++it) {
                const libMesh::Elem *ele = *it;

                for (int kk = 0; kk < ele->n_sides(); kk++) {
                    auto neigh = ele->neighbor_ptr(kk);

                    auto bc_id = utopia::boundary_id(mesh.get_boundary_info(), ele, kk);

                    auto check = (std::find(dirichlet_id.begin(), dirichlet_id.end(), bc_id) != dirichlet_id.end());

                    if (check) {
                        index_local.clear();

                        auto side = ele->build_side_ptr(kk);

                        for (int ll = 0; ll < side->n_nodes(); ll++) {
                            const libMesh::Node *node = side->node_ptr(ll);

                            const libMesh::dof_id_type node_dof = node->dof_number(sys_number, var_number, 0);

                            if (on_boundary.count(node->id()) && dof_map.is_constrained_dof(node_dof))
                                index.push_back(node_dof);
                        }
                    }
                }
            }
        }

     //   std::cout << "Adaptivity::compute_boundary_nodes::END " << std::endl;
       #endif
       std::cout << "Adaptivity::compute_boundary_nodes::END " << std::endl;
    }

    void Adaptivity::compute_boundary_nodes_to_skip(const libMesh::MeshBase &mesh,
                                                    libMesh::DofMap &dof_map,
                                                    unsigned int sys_number,
                                                    unsigned int var_number,
                                                    std::vector<int> & index)
   {
     std::cout << "Adaptivity::compute_boundary_nodes_to_skip::BEGIN " << std::endl;
     
     #ifdef LIBMESH_VERSION_LESS_THAN(1, 4, 0)
     
     #else
        // Only constrain elements in 2,3D.
        if (mesh.mesh_dimension() == 1) return;

        unsigned int mesh_dim = mesh.mesh_dimension();

        // std::cout<<"lagrange_compute_constraints libmesh mio
        // prima:"<<constraints.size()<<std::endl;
        // libmesh_assert(elem);

        std::vector<int> index_self;

        index_self.clear();

        unsigned int var_num = 0;

        libMesh::MeshBase::const_element_iterator it_0 = mesh.active_elements_begin();

        const libMesh::MeshBase::const_element_iterator end_it_0 = mesh.active_elements_end();

        for (; it_0 != end_it_0; ++it_0) {
            const libMesh::Elem *elem = *it_0;

            // Only constrain active and ancestor elements
            if (elem->subactive())  // if the element is subactive (i.e. has no active
                                    // descendants)
                return;

            libMesh::FEType fe_type = dof_map.variable_type(var_num);
            fe_type.order = static_cast<libMesh::Order>(fe_type.order + elem->p_level());

            std::vector<libMesh::dof_id_type> my_dof_indices, parent_dof_indices;
            libMesh::UniquePtr<const libMesh::Elem> my_side, parent_side;

            for (auto s : elem->side_index_range()) {
                if (elem->neighbor_ptr(s) != nullptr && elem->neighbor_ptr(s) != libMesh::remote_elem) {
                    if (elem->neighbor_ptr(s)->level() < elem->level())

                    {
                        const auto *parent = elem->parent();
                        libmesh_assert(parent);

                        elem->build_side_ptr(my_side, s);

                        parent->build_side_ptr(parent_side, s);

                        my_dof_indices.reserve(my_side->n_nodes());

                        parent_dof_indices.reserve(parent_side->n_nodes());

                        dof_map.dof_indices(my_side.get(), my_dof_indices, var_num);

                        dof_map.dof_indices(parent_side.get(), parent_dof_indices, var_num);

                        const unsigned int n_side_dofs =
                            libMesh::FEInterface::n_dofs(mesh_dim - 1, fe_type, my_side->type());

                        const unsigned int n_parent_side_dofs =
                            libMesh::FEInterface::n_dofs(mesh_dim - 1, fe_type, parent_side->type());

                        for (unsigned int my_dof = 0; my_dof != n_side_dofs; my_dof++) {
                            libmesh_assert_less(my_dof, my_side->n_nodes());
                            assert(my_dof < n_side_dofs);

                            const libMesh::dof_id_type my_dof_g = my_dof_indices[my_dof];

                            bool self_constraint = false;

                            for (unsigned int their_dof = 0; their_dof != n_parent_side_dofs; their_dof++) {
                                libmesh_assert_less(their_dof, parent_side->n_nodes());

                                const libMesh::dof_id_type their_dof_g = parent_dof_indices[their_dof];

                                if (their_dof_g == my_dof_g) {
                                    // self_constraint = true;

                                    auto check_2 = (std::find(index_self.begin(), index_self.end(), their_dof_g) !=
                                                    index_self.end());

                                    if (!check_2) {
                                        index_self.push_back(their_dof_g);

                                        // std::cout<<"their_dof_g"<<their_dof_g<<std::endl;
                                    }

                                    // index_self.push_back(their_dof_g);

                                    // std::cout<<"their_dof_g"<<their_dof_g<<std::endl;
                                    // break;
                                }
                            }
                        }
                    }
                }
            }
        }
        // std::cout<<"Adaptivity::compute_boundary_nodes::Begin "<<std::endl;

        auto on_boundary = libMesh::MeshTools::find_boundary_nodes(mesh);

        std::vector<int> dirichlet_id, index_local, tmp;

        index_local.clear();

        dirichlet_id.clear();

        index.clear();

        libMesh::MeshBase::const_element_iterator it_1 = mesh.active_elements_begin();

        const libMesh::MeshBase::const_element_iterator end_it_1 = mesh.active_elements_end();

        libMesh::UniquePtr<const libMesh::Elem> parent_side_0_new;

        for (; it_1 != end_it_1; ++it_1) {
            const libMesh::Elem *ele = *it_1;

            const auto *ele_parent_0 = ele->top_parent();

            for (int jj = 0; jj < ele_parent_0->n_sides(); jj++) {
                libmesh_assert(ele_parent_0);

                auto parent_side_0_new = ele_parent_0->build_side_ptr(jj);

                index_local.clear();

                for (int ll = 0; ll < parent_side_0_new->n_nodes(); ll++) {
                    const libMesh::Node *node_0 = parent_side_0_new->node_ptr(ll);

                    const libMesh::dof_id_type node_dof_0 = node_0->dof_number(sys_number, var_number, 0);

                    if (dof_map.is_constrained_dof(node_dof_0)) {
                        index_local.push_back(node_dof_0);
                    }
                }

                if (index_local.size() == parent_side_0_new->n_nodes()) {
                    auto bc_id = utopia::boundary_id(mesh.get_boundary_info(), ele_parent_0, jj);

                    auto check = (std::find(dirichlet_id.begin(), dirichlet_id.end(), bc_id) != dirichlet_id.end());

                    if (!check) dirichlet_id.push_back(bc_id);
                }
            }
        }

        libMesh::MeshBase::const_element_iterator it_2 = mesh.active_elements_begin();

        const libMesh::MeshBase::const_element_iterator end_it_2 = mesh.active_elements_end();

        for (; it_2 != end_it_2; ++it_2) {
            const libMesh::Elem *ele = *it_2;

            // if (ele->subactive())
            {  // if the element is subactive (i.e. has no active descendants)

                for (int kk = 0; kk < ele->n_sides(); kk++) {
                    auto neigh = ele->neighbor_ptr(kk);
                    if (ele->neighbor_ptr(kk) == nullptr && ele->neighbor_ptr(kk) != libMesh::remote_elem) {
                        auto bc_id = utopia::boundary_id(mesh.get_boundary_info(), ele, kk);

                        auto check = (std::find(dirichlet_id.begin(), dirichlet_id.end(), bc_id) != dirichlet_id.end());

                        if (check) {
                            index_local.clear();

                            auto side = ele->build_side_ptr(kk);

                            for (int ll = 0; ll < side->n_nodes(); ll++) {
                                const libMesh::Node *node = side->node_ptr(ll);

                                const libMesh::dof_id_type node_dof = node->dof_number(sys_number, var_number, 0);

                                if (on_boundary.count(node->id()) && dof_map.is_constrained_dof(node_dof)) {
                                    auto check_2 = (std::find(tmp.begin(), tmp.end(), node_dof) != tmp.end());

                                    if (!check_2) {
                                        tmp.push_back(node_dof);
                                    }
                                    // std::cout<<"tmp"<<node_dof<<std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }

        for (auto it = tmp.begin(); it != tmp.end(); ++it) {
            auto check = (std::find(index_self.begin(), index_self.end(), *it) != index_self.end());

            if (!check) {
                auto check_2 = (std::find(index.begin(), index.end(), *it) != index.end());

                if (!check_2) {
                    index.push_back(*it);
                    // std::cout<<"index_to_skip"<<*it<<std::endl;
                }
            }
        }

        std::cout << "Adaptivity::compute_boundary_nodes_to_skip::END " << tmp.size() << " and " << index.size()
                 << std::endl;
    #endif
    }
}  // namespace utopia

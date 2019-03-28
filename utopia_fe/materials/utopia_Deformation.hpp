#ifndef UTOPIA_DEFORMATION_HPP
#define UTOPIA_DEFORMATION_HPP

#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"

#include "utopia.hpp"
#include "utopia_fe_base.hpp"

#include <vector>
#include <algorithm>



namespace utopia {

    inline static void deform_mesh(
        libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map,
        const UVector &disp,
        const std::vector<int> &vars_in = std::vector<int>(),
        int sys_num = 0)
    {

        const auto dim = mesh.mesh_dimension();

        std::vector<int> vars;
        if(vars_in.empty()) {
            vars.resize(dim);
            int n = 0;
            std::generate(vars.begin(), vars.end(), [&n] () { return n++; });
        } else {
            vars = vars_in;
        }



        if(mesh.comm().size() > 1) {
            auto r = range(disp);
            Read<UVector> r_d(disp);

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
                        idx_to_value[dof_id] = disp.get(dof_id);
                    } else {
                        unique_idx.insert(dof_id);
                    }
                }
            }

            idx.insert(idx.end(), unique_idx.begin(), unique_idx.end());
            UVector out = disp.select(idx);
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

                    for(unsigned int c = 0; c < vars.size(); ++c) {
                        const int dof_id = node.dof_number(sys_num, vars[c], 0);
                        assert(idx_to_value.find(dof_id) != idx_to_value.end());
                        double &val = idx_to_value[dof_id];
                        node(c) += val;
                        val = 0.;
                    }
                }
            }

        } else {

            Read<UVector> r_d(disp);

            auto m_it  = mesh.local_nodes_begin();
            auto m_end = mesh.local_nodes_end();


            for(; m_it != m_end; ++m_it) {
                for(unsigned int c = 0; c < vars.size(); ++c) {
                    const int dof_id = (*m_it)->dof_number(sys_num, vars[c], 0);
                    (**m_it)(c) += disp.get(dof_id);
                }
            }
        }
    }

}

#endif //UTOPIA_DEFORMATION_HPP

#include "utopia_NormalTangentialCoordinateSystem.hpp"

#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/libmesh.h"
#include "libmesh/quadrature_gauss.h"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_intersector.hpp"

// #include "MortarAssemble.hpp"
// #include "utopia_Socket.hpp"

#include <cmath>

namespace utopia {
    void scale_normal_vector_with_gap(const int dim, const DVectord &normals, const DVectord &gap, DVectord &out)
    {
        out = local_zeros(local_size(normals));

        Read<DVectord> r_n(normals), r_g(gap);
        Write<DVectord> w_o(out);

        using namespace libMesh;

        Range r = range(normals);

        std::vector<double> n(dim);
        for(auto i = r.begin(); i < r.end(); i += dim) {

            double len = 0.;
            for(int d = 0; d < dim; ++d) {
                n[d] = normals.get(i + d);
                len += n[d] * n[d];
            }

            len = std::sqrt(len);

            if(len < 1e-8) {
                continue;
            }

            double g = gap.get(i);
            for(int d = 0; d < dim; ++d) {
                out.set(i+d, n[d] * g);
            }
        }
    }



    
    bool assemble_normal_tangential_transformation(const libMesh::MeshBase &mesh,
                                                   const libMesh::DofMap &dof_map,
                                                   const std::vector<int> &boundary_tags,
                                                   DVectord &is_normal_component,
                                                   DVectord &normals,
                                                   DSMatrixd &mat)
    {
        using namespace libMesh;
        
        Intersector isector;
        
        const uint n_dims = mesh.mesh_dimension();
        std::unique_ptr<FEBase> fe = FEBase::build(n_dims, dof_map.variable_order(0));
        fe->get_normals();
        fe->get_phi();
        fe->get_JxW();
        
        QGauss quad(n_dims-1 , FIFTH);
        
        
        const SizeType local_dofs = dof_map.n_local_dofs();
        DVectord global_normal_vec = local_zeros(local_dofs);
        DVectord touched = local_zeros(local_dofs);

        normals = local_zeros(local_dofs);
        
        std::vector<libMesh::dof_id_type> dof_indices;
        DenseVector<Real> vec, local_touched;
        
        SizeType n_detected_side_sets = 0;
        { //synch-block begin
            Write<DVectord> w_n(global_normal_vec), w_t(touched);
            
            for(auto e_it = mesh.active_local_elements_begin();
                e_it != mesh.active_local_elements_end(); ++e_it) {
                const auto &e = **e_it;

                for(uint side = 0; side < e.n_sides(); ++side) {
                    if(e.neighbor_ptr(side) != nullptr) {continue;}
                    
                    bool select = false;
                    for(auto t : boundary_tags) {
                        if(mesh.get_boundary_info().has_boundary_id(&e, side, t)) {
                            select = true;
                            break;
                        }
                    }
                    
                    if(!select) continue;
                    ++n_detected_side_sets;
                    
                    fe->attach_quadrature_rule(&quad);
                    fe->reinit(&e, side);
                    
                    const auto &fe_normals = fe->get_normals();
                    const auto &fun        = fe->get_phi();
                    const auto &JxW        = fe->get_JxW();
                    
                    const uint n_fun = fun.size();
                    const uint n_qp  = fun[0].size();
                    
                    vec.resize(n_fun * n_dims);
                    vec.zero();
                    
                    for(uint qp = 0; qp < quad.n_points(); ++qp) {
                        for(uint i = 0; i < fun.size(); ++i){
                            for(uint d = 0; d < n_dims; ++d) {
                                vec(i + d * n_fun) += fun[i][qp] * fe_normals[qp](d) * JxW[qp];
                            }
                        }
                    }
                    
                    dof_map.dof_indices(&e, dof_indices);
                    
                    assert(dof_indices.size() == n_fun * n_dims);
                    
                    for(uint i = 0; i < dof_indices.size(); ++i) {
                        const uint ind = dof_indices[i];
                        global_normal_vec.add(ind, vec(i));
                        touched.set(ind, std::abs(vec(i)) > 1e-16);
                    }
                }
            }
        } //synch-block end
        
        // write("v.m", global_normal_vec);
        mat = local_sparse(local_dofs, local_dofs, n_dims);        
        auto r = range(global_normal_vec);
        
        std::vector<Real> H(n_dims*n_dims, 0);
        
        is_normal_component = local_zeros(local_dofs);
        
        SizeType n_detecetd_normals = 0;
        { //synch-block begin
            Read<DVectord> r_n(global_normal_vec), r_t(touched);
            Write<DSMatrixd> w_m(mat);
            Write<DVectord> w_i(is_normal_component);
            Write<DVectord> w_n(normals);
            
            for(SizeType i = r.begin(); i < r.end(); i += n_dims) {
                bool is_node_touched = false;
                for(uint j = 0; j < n_dims; ++j) {
                   is_node_touched = is_node_touched || touched.get(i+j) > 0;
                }
                
                bool use_identity = false;
                
                if(!is_node_touched) {
                    use_identity = true;
                } else {
                    ++n_detecetd_normals;

                    std::vector<Real> n(n_dims, 0);
                    Real norm = 0.;

                    for(uint j = 0; j < n_dims; ++j) {
                        n[j] = global_normal_vec.get(i + j);
                        norm += n[j] * n[j];
                    }


                    is_normal_component.set(i, 1.);
                    
                    norm = std::sqrt(norm);
                    
                    for(uint j = 0; j < n_dims; ++j) {
                        n[j] /= norm;

                        normals.set(i + j, n[j]);
                    }


                    
                    n[0] -= 1.;
                    
                    if(std::abs(n[0]) < 1e-16) {
                        use_identity = true;
                    } else {
                        use_identity = false;

                        if(n_dims == 2) {
                            norm = std::sqrt(n[0] * n[0] + n[1] * n[1]);
                            n[0] /= norm; n[1] /= norm;

                            isector.householder_reflection_2(&n[0], &H[0]);
                        } else {
                            norm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                            n[0] /= norm; n[1] /= norm; n[2] /= norm;

                            isector.householder_reflection_3(&n[0], &H[0]);
                        }
                        
                        for(uint di = 0; di < n_dims; ++di) {
                            for(uint dj = 0; dj < n_dims; ++dj) {
                                mat.set((i + di), (i + dj), H[di * n_dims + dj]);
                            }
                        }
                    }
                }
                
                if(use_identity) {
                    for(uint di = 0; di < n_dims; ++di) {
                        mat.set((i + di), (i + di), 1.0);
                    }
                }
            }
        } //synch-block end
        
        return true;
    }
    
}

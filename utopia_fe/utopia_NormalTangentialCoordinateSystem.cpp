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
#include "utopia_fe.hpp"
#include "utopia_LibMeshBackend.hpp"

#include "MortarAssemble.hpp"
#include "utopia_Socket.hpp"

#include <cmath>

namespace utopia {


	void plot_scaled_normal_field(const libMesh::MeshBase &mesh,
								  const DVectord &normals,
								  const DVectord &scale,
								  const std::string &name)
	{
#ifdef WITH_BOOST
		using namespace libMesh;
		int mesh_dim = mesh.mesh_dimension();
		

		DenseVector<double> local_normal;
		DenseVector<Real> local_scale;
		
		std::vector<double> all_points, all_normals;
		
		std::vector<double> point(mesh_dim, 0.);
		for(auto n_it = mesh.active_nodes_begin(); n_it != mesh.active_nodes_end(); ++n_it) {
			Node &n = **n_it;
			
			std::vector<dof_id_type> node_dof_ids;
			
			for(int d = 0; d < mesh_dim; ++d) {
				auto dof_id = n.dof_number(0, d, 0);
				node_dof_ids.push_back(dof_id);
				
				point[d] = n(d);
			}
			
			get_vector(normals, node_dof_ids, local_normal);
			get_vector(scale,   node_dof_ids, local_scale);
			
			for(int d = 0; d < mesh_dim; ++d) {
				local_normal(d) *= local_scale(0);
			}

			if(local_normal.l2_norm() < 1e-16) continue;
			
			all_points.insert(all_points.end(),
							  point.begin(),
							  point.end());

			all_normals.insert(all_normals.end(),
							   local_normal.get_values().begin(),
							   local_normal.get_values().end());
		}
		
		if(all_points.empty()) {
			return;
		}
		
		quiver(mesh_dim,
			   all_points.size()/mesh_dim,
			   &all_points[0],
			   &all_normals[0],
			   name);
#endif //WITH_BOOST
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
        std::unique_ptr<FEBase> fe = FEBase::build(n_dims, FIRST);
        
        QGauss quad(n_dims-1 , FIFTH);
        
        
        const SizeType local_dofs = dof_map.n_local_dofs();
        DVectord global_normal_vec = local_zeros(local_dofs);
        normals = local_zeros(local_dofs);
        
        std::vector<libMesh::dof_id_type> dof_indices;
        DenseVector<Real> vec;
        
        SizeType n_detected_side_sets = 0;
        { //synch-block begin
            Write<DVectord> w_n(global_normal_vec);
            
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
                    
                    //assemble weighted normal
                    const auto &fun = fe->get_phi();
                    const auto &JxW = fe->get_JxW();
                    
                    const uint n_fun = fun.size();
                    const uint n_qp  = fun[0].size();
                    
                    vec.resize(n_fun * n_dims);
                    vec.zero();
                    
                    for(uint qp = 0; qp < quad.n_points(); ++qp) {
                        for(uint i = 0; i < fun.size(); ++i){
                            for(uint d = 0; d < n_dims; ++d) {
                                vec(i * n_dims + d) += fun[i][qp] * fe_normals[i](d) * JxW[qp];
                            }
                        }
                    }
                    
                    dof_map.dof_indices(&e, dof_indices);
                    
                    assert(dof_indices.size() == n_fun * n_dims);
                    
                    for(uint i = 0; i < dof_indices.size(); ++i) {
                        const uint ind = dof_indices[i];
                        global_normal_vec.add(ind, vec(i));
                    }
                }
            }
        } //synch-block end

        disp(global_normal_vec);
        
        mat = local_sparse(local_dofs, local_dofs, n_dims);
        
        auto r = range(global_normal_vec);
        
        std::vector<Real> H(n_dims*n_dims, 0);
        
        is_normal_component = local_zeros(local_dofs * n_dims);
        
        SizeType n_detecetd_normals = 0;
        { //synch-block begin
            Read<DVectord> r_n(global_normal_vec);
            Write<DSMatrixd> w_m(mat);
            Write<DVectord> w_i(is_normal_component);
            Write<DVectord> w_n(normals);
            
            for(SizeType i = r.begin(); i < r.end(); i += n_dims) {
                std::vector<Real> n(n_dims, 0);
                Real norm = 0.;

                for(uint j = 0; j < n_dims; ++j) {
                    n[j] = global_normal_vec.get(i + j);
                    norm += n[j] * n[j];
                }
                
                bool use_identity = false;
                
                if(norm < 1e-16) {
                    use_identity = true;
                } else {
                    ++n_detecetd_normals;
                    is_normal_component.set(i, 1.);
                    
                    norm = std::sqrt(norm);
                    
                    for(uint j = 0; j < n_dims; ++j) {
                        n[j] /= norm;

                        normals.set(i + j, n[j]);
                    }


                    
                    n[0] -= 1;
                    
                    if(std::abs(n[0]) < 1e-16) {
                        use_identity = true;
                    } else {
                        if(n_dims == 2) {
                            isector.householder_reflection_2(&n[0], &H[0]);
                        } else {
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
        
        // std::cout << "n_detected_side_sets: " << n_detected_side_sets << std::endl;	
        // std::cout << "n_detecetd_normals: " << n_detecetd_normals << std::endl;


        // plot_scaled_normal_field(mesh,
								//  normals,
								//   local_values(local_dofs, 0.1),
								//   "normals");
        
       
        return true;
    }
    
}

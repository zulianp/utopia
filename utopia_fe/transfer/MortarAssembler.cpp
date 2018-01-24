#include "MortarAssembler.hpp"
#include "MortarAssemble.hpp"
#include "utopia_triangulate.hpp"

#include <cmath>
#include <queue>
#include <numeric>

namespace utopia {

	MortarAssembler::MortarAssembler(
		const std::shared_ptr<LibMeshFESpaceBase> &master, 
		const std::shared_ptr<LibMeshFESpaceBase> &slave)
	: master_(master), slave_(slave), use_biorth_(false)
	{ }


	bool MortarAssembler::assemble(DSMatrixd &B)
	{
		using namespace std;
		static const bool verbose = true;

		const auto &master_mesh = master_->mesh();
		const auto &slave_mesh  = slave_->mesh();

		int dim = master_mesh.mesh_dimension();

		// Chrono c;
		// c.start();

		std::vector<int> pairs;
		if(!hash_grid_detect_intersections(master_mesh, slave_mesh, pairs)) {
			return false;
		}

		// c.stop();

		if(verbose) {
			std::cout << "intersection detection:\t";
			// c.describe(std::cout);
		}


		// c.start();

		libMesh::DenseMatrix<libMesh::Real> master_pts;
		libMesh::DenseMatrix<libMesh::Real> slave_pts;
		libMesh::DenseMatrix<libMesh::Real> intersection2;

		Polyhedron master_poly, slave_poly;
		Polyhedron  intersection3;
		Intersector isector;

		std::shared_ptr<Transform> master_trans;
		std::shared_ptr<Transform> slave_trans;

		std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;

		master_fe = libMesh::FEBase::build(master_mesh.mesh_dimension(), 
			master_->dof_map().variable_type(master_->var_num()));

		slave_fe  = libMesh::FEBase::build(slave_mesh.mesh_dimension(), 
			slave_->dof_map().variable_type(slave_->var_num()));

		const int master_order = master_->order();
		const int slave_order  = slave_->order();


	 	//////////////////////////////////////////////////
//		int skip_zeros = 1;
		B = sparse(slave_->dof_map().n_local_dofs(),
			master_->dof_map().n_local_dofs(), 
			std::max(1, int(master_->dof_map().n_local_dofs() * 0.2)));

		D = sparse(slave_->dof_map().n_local_dofs(), slave_->dof_map().n_local_dofs(), std::max(1, int(master_->dof_map().n_local_dofs() * 0.2)));
		
		{ //write scope begin
			Write<DSMatrixd> w_B(B);
			Write<DSMatrixd> w_D(D);

			std::vector<libMesh::dof_id_type> master_dofs, slave_dofs;
			libMesh::DenseMatrix<libMesh::Real> elemmat;
			libMesh::DenseMatrix<libMesh::Real> other_mat;
			libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;

			QMortar composite_ir(dim);
			QMortar master_ir(dim);
			QMortar slave_ir(dim);
	 	// //////////////////////////////////////////////////

			libMesh::Real total_intersection_volume = 0.0;
			libMesh::Real local_element_matrices_sum = 0.0;
			libMesh::Real other_sum = 0.0;

			std::vector<libMesh::Real> volumes(slave_mesh.n_active_local_elem());

			int i = 0;
			for(auto e_it = slave_mesh.active_local_elements_begin(); e_it != slave_mesh.active_local_elements_end(); ++e_it, ++i) {
				const auto &elem = **e_it;
				if(dim == 2) {
					make_polygon(elem, slave_pts);
					volumes[i] = isector.polygon_area_2(slave_pts.m(), &slave_pts.get_values()[0]);
				} else if(dim == 3) {
					make_polyhedron(elem, intersection3);
					volumes[i] = isector.p_mesh_volume_3(intersection3);
				} else {
					assert(false);
				}
			}

			bool intersected = false;
			long n_intersections = 0;
			for(auto it = begin(pairs); it != end(pairs);    ) {
				const int master_index = *it++;
				const int slave_index  = *it++;

				auto &master_el = *master_mesh.elem(master_index);
				auto &slave_el  = *slave_mesh.elem(slave_index);

				libMesh::Real weight   = volumes[slave_index];

				const int order = order_for_l2_integral(dim, master_el, master_order, slave_el, slave_order);

				bool pair_intersected = false;
				if(dim == 2)  {	
					make_polygon(master_el, master_pts);
					make_polygon(slave_el, slave_pts);

					if(slave_el.has_affine_map() && master_el.has_affine_map()) {
						if(intersect_2D(master_pts, slave_pts, intersection2)) {
							const double area = isector.polygon_area_2(intersection2.m(), &intersection2.get_values()[0]);
							total_intersection_volume += fabs(area);
							make_composite_quadrature_2D(intersection2, weight, order, composite_ir);

							pair_intersected = true;

							master_trans = std::make_shared<Transform2>(master_el);
							slave_trans  = std::make_shared<Transform2>(slave_el);
						} 
					} else {
						std::vector<int>    slave_tri, master_tri;
						std::vector<int>    intersection_tri;
						std::vector<double> intersection_points;
						std::vector<int> temp_tri;

						triangulate_polygon(slave_pts.m(), &slave_pts.get_values()[0],   slave_tri);
						triangulate_polygon(master_pts.m(), &master_pts.get_values()[0],  master_tri);

						double relative_intersection_volume = 0;
						libMesh::DenseMatrix<libMesh::Real> m_tri(3, 2), s_tri(3, 2);
						for(size_t i_slave = 0; i_slave < slave_tri.size(); i_slave += 3) {
							const int v1_s = slave_tri[i_slave];
							const int v2_s = slave_tri[i_slave + 1];
							const int v3_s = slave_tri[i_slave + 2];

							s_tri(0, 0) = slave_pts(v1_s, 0); s_tri(0, 1) = slave_pts(v1_s, 1); 
							s_tri(1, 0) = slave_pts(v2_s, 0); s_tri(1, 1) = slave_pts(v2_s, 1); 
							s_tri(2, 0) = slave_pts(v3_s, 0); s_tri(2, 1) = slave_pts(v3_s, 1); 

							for(size_t i_master = 0; i_master < master_tri.size(); i_master += 3) {
								const int v1_m = master_tri[i_master];
								const int v2_m = master_tri[i_master + 1];
								const int v3_m = master_tri[i_master + 2];

								m_tri(0, 0) = master_pts(v1_m, 0); m_tri(0, 1) = master_pts(v1_m, 1); 
								m_tri(1, 0) = master_pts(v2_m, 0); m_tri(1, 1) = master_pts(v2_m, 1); 
								m_tri(2, 0) = master_pts(v3_m, 0); m_tri(2, 1) = master_pts(v3_m, 1); 

								assert(isector.polygon_area_2(m_tri.m(), &m_tri.get_values()[0]) > 0);
								assert(isector.polygon_area_2(s_tri.m(), &s_tri.get_values()[0]) > 0);

								if(intersect_2D(m_tri, s_tri, intersection2)) {
									const double area = isector.polygon_area_2(intersection2.m(), &intersection2.get_values()[0]);
									relative_intersection_volume += fabs(area);

									triangulate_polygon(intersection2.m(), &intersection2.get_values()[0],  temp_tri);

									const int offset = intersection_points.size()/2;
									std::transform(temp_tri.begin(), temp_tri.end(), temp_tri.begin(), 
										[offset](const int val) -> int { return offset + val; });

									intersection_tri.insert(intersection_tri.end(), temp_tri.begin(), temp_tri.end());
									intersection_points.insert(intersection_points.end(), intersection2.get_values().begin(), intersection2.get_values().end());
									pair_intersected = true;
								}
							}
						}

						if(pair_intersected) {
							total_intersection_volume += relative_intersection_volume;
							libMesh::DenseMatrix<libMesh::Real> isect_mesh_pts(intersection_points.size()/2, 2);
							isect_mesh_pts.get_values() = std::move(intersection_points);
							make_composite_quadrature_2D_from_tri_mesh(intersection_tri, isect_mesh_pts, weight, order, composite_ir);

							master_trans = std::make_shared<Transform2>(master_el);
							slave_trans  = std::make_shared<Transform2>(slave_el);
						}
					}

				} else if(dim == 3) {
					make_polyhedron(master_el, master_poly);
					make_polyhedron(slave_el,  slave_poly);

					if(intersect_3D(master_poly, slave_poly, intersection3)) {
						total_intersection_volume += isector.p_mesh_volume_3(intersection3);

	 				// const int order = master_order * ( is_hex(master_el.type())? 3 : 1 ) + 
	 				// 				  slave_order  * ( is_hex(slave_el.type())? 3  : 1 ) * (slave_el.has_affine_map()? 1 : 2);

						make_composite_quadrature_3D(intersection3, weight, order, composite_ir);
						pair_intersected = true;

						master_trans = std::make_shared<Transform3>(master_el);
						slave_trans  = std::make_shared<Transform3>(slave_el);
					}

				} else {
					assert(false);
					return false;
				}

				if(pair_intersected) {
	 				//make reference quaratures
					transform_to_reference(*master_trans, master_el.type(), composite_ir,  master_ir);
					transform_to_reference(*slave_trans,  slave_el.type(),  composite_ir,  slave_ir);

					master_->dof_map().dof_indices(&master_el, master_dofs);
					slave_->dof_map().dof_indices(&slave_el,   slave_dofs);

					master_fe->attach_quadrature_rule(&master_ir);
					master_fe->reinit(&master_el);

					slave_fe->attach_quadrature_rule(&slave_ir);
					slave_fe->reinit(&slave_el);

					elemmat.zero();
					other_mat.zero();

					if(use_biorth_) {
						// std::unique_ptr<libMesh::FEBase> biorth_elem = 
						// libMesh::FEBase::build(slave_->mesh().mesh_dimension(), 
						// 					   slave_->dof_map().variable_type(slave_->var_num()));

						// libMesh::QGauss qg(dim, libMesh::Order(order));
						// biorth_elem->attach_quadrature_rule(&qg);
						// biorth_elem->reinit(&slave_el);

						// libMesh::DenseMatrix<libMesh::Real> weights;
						// mortar_assemble_weights(*biorth_elem, weights);

						mortar_assemble_biorth(*master_fe, *slave_fe, slave_el.type(), elemmat);
						mortar_assemble_biorth(*slave_fe,  *slave_fe, slave_el.type(), other_mat);

					} else {
						mortar_assemble(*master_fe, *slave_fe, elemmat);
						mortar_assemble(*slave_fe, *slave_fe, other_mat);
					}
					

					add_matrix(elemmat,   slave_dofs, master_dofs, B);
					add_matrix(other_mat, slave_dofs, slave_dofs, D);

					local_element_matrices_sum += std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(), libMesh::Real(0.0));
					other_sum += std::accumulate(other_mat.get_values().begin(), other_mat.get_values().end(), libMesh::Real(0.0));
					intersected = true;
					++n_intersections;
				}
			} 

			if(!intersected) return false;

			// c.stop();

			if(verbose) {	
				std::cout << "assembly time:\t";
				// c.describe(std::cout);

				std::cout << "n_intersections:   " << n_intersections << std::endl;
				std::cout << "n_false_positives: " << (pairs.size()/2 - n_intersections) << std::endl;
				std::cout << "local_element_matrices_sum: " << local_element_matrices_sum << std::endl;
				std::cout << "intersection volume: " << total_intersection_volume << std::endl;
				std::cout << "B in R^(" << size(B).get(0) <<  " x " << size(B).get(1) << ")" << std::endl;

				std::cout << "other_sum: " << other_sum << std::endl;
			}

		} //write scope end

		// write("t_D.m", D);
		// write("t_B.m", B);

		return true;
	}


	void ContactAssembly::finalize()
	{
		//Compute average gap
		const auto &v = gap.get_values();
		const auto &m = coupling.get_values();
		avg_gap = (std::accumulate(v.begin(), v.end(), RealT(0.))/std::accumulate(m.begin(), m.end(), RealT(0.)))*normals.n();
	}

	MortarContactAssembler::MortarContactAssembler(const std::shared_ptr<LibMeshFESpaceBase> &space)
	: space_(space), strict_gap_policy(false)
	{}

	void build_adj_list(const libMesh::MeshBase &mesh, std::vector< std::vector<long> > &adj_list)
	{
		std::vector< std::vector<long> > vertex_to_boundary_element(mesh.n_local_nodes());

		uint max_id = 0;

		for(auto e_it = mesh.active_local_elements_begin(); e_it != mesh.active_local_elements_end(); ++e_it) {
			auto &el_ptr = *e_it;

			for(uint side = 0; side < el_ptr->n_sides(); ++side) {
				if(el_ptr->neighbor_ptr(side) != nullptr) continue;
				
				auto side_ptr = el_ptr->build_side_ptr(side);	

				for(uint node = 0; node < side_ptr->n_nodes(); ++node) {
					vertex_to_boundary_element[side_ptr->node_id(node)].push_back(el_ptr->id());
					max_id = std::max(max_id, el_ptr->id());
				}
			}
		}

		adj_list.resize(max_id + 1);

		for(const auto &v2f : vertex_to_boundary_element) {
			for(uint i = 0; i < v2f.size(); ++i) {
				for(uint j = 0; j < v2f.size(); ++j) {
					if(v2f[i] != v2f[j]) {
						adj_list[v2f[i]].push_back(v2f[j]);
					}
				}
			}
		}

		for(auto &adj : adj_list) {
			std::sort(adj.begin(),   adj.end());
			auto it = std::unique(adj.begin(), adj.end());
			adj.erase(it, adj.end());

			// std::cout << adj.size() << " : " << std::endl;
			// print_vector(adj.begin(), adj.end());
			// std::cout << "-----------------------\n";
		}
	}

	void build_dag(std::vector< std::shared_ptr<ContactAssembly> > &contacts, std::vector< std::vector<long> > &dag, std::vector<long> &ordering)
	{
		long max_side_id = 0;
		for(auto c_ptr : contacts) {
			max_side_id = std::max(c_ptr->id_slave,  max_side_id);
			// c_ptr->describe(std::cout);
		}

		std::vector<libMesh::Real> areas(max_side_id + 1, 0);
		dag.resize(areas.size());

		for(auto c_ptr : contacts) {
			areas[c_ptr->id_slave] +=  c_ptr->relative_area;
			c_ptr->is_valid = true;
		}

		for(auto c_ptr : contacts) {
			if(std::abs(areas[c_ptr->id_slave]) + 1e-8 < 1) {
				//the element is eliminated as a slave
				c_ptr->is_valid = false;
			}
		}

		std::sort(contacts.begin(), contacts.end(), [](const std::shared_ptr<ContactAssembly> &left, const std::shared_ptr<ContactAssembly> &right) -> bool {
			return *left < *right;
		});

		ordering.clear();

		std::vector<bool> visited(max_side_id + 1, false);

		for(auto c_ptr : contacts) {
			if(c_ptr->is_valid) {
				dag[c_ptr->id_slave].push_back(c_ptr->id_master);
				if(!visited[c_ptr->id_slave]) {
					ordering.push_back(c_ptr->id_slave);
					visited[c_ptr->id_slave] = true;
				}
			}
		}
	}

	void breadth_first_color(const long seed, 
		const ushort color,
		const std::vector< std::vector<long> > &adj_list, 
		const std::vector< std::vector<long> > &dag, 
		std::vector<ushort> &role)
	{
		std::queue<long> boundary_elements;
		boundary_elements.push(seed);

		role[seed] = color;

		while (!boundary_elements.empty()) {
			const long current = boundary_elements.front();
			boundary_elements.pop();
			const auto &adj = adj_list[current];

			for (uint i = 0; i < adj.size(); ++i) {
				const long next = adj[i];

				if (next < dag.size() && !dag[next].empty()) {
					if(role[next] == UNASSIGNED) {
						role[next] = color;
						boundary_elements.push(next);
					} else {
						// assert(role[next] == color);
					}
				}
			}
		}
	}

	void assign_master_and_slave_roles(const std::vector< std::vector<long> > &dag, const std::vector<long> &ordering, const std::vector< std::vector<long> > &adj_list, std::vector<ushort> &role)
	{
		size_t source = 0;

		std::vector<size_t> queue;
		queue.reserve(20);

		role.resize(dag.size()); 
		std::fill(role.begin(), role.end(), UNASSIGNED);

		for(uint i = 0; i < ordering.size(); ++i) {
			const long el = ordering[i];

			// if(role[el] != UNASSIGNED) continue;
			if(dag[el].empty()) continue;

			queue.push_back(el);

			while(!queue.empty()) {
				source = queue.back();
				queue.pop_back();

				if(role[source] == UNASSIGNED) {
					role[source] = SLAVE;
				}

				bool must_change_role_to_removed = false;
				const auto &paired = dag[source];

				if(role[source] == SLAVE) {
					for(auto a : paired) {
						if(a < role.size() && role[a] == SLAVE) {
							must_change_role_to_removed = true;
							break;
						}
					}

					if(must_change_role_to_removed) {
						std::cout << "WARNING SLAVE CHANGING TO REMOVED\n" << std::endl;
						role[source] = REMOVED;

					} else {
						breadth_first_color(source, SLAVE, adj_list, dag, role);

						for(auto a : paired) {
							if(a < role.size() && role[a] == UNASSIGNED) {
								queue.push_back(a);
								role[a] = MASTER;

								breadth_first_color(a, MASTER, adj_list, dag, role);
							}
						}
					}
				}  

				if(role[source] == MASTER) {
					breadth_first_color(source, MASTER, adj_list, dag, role);

					for(auto a : paired) {
						if(a < role.size() && role[a] == UNASSIGNED) {
							queue.push_back(a);

							breadth_first_color(a, SLAVE, adj_list, dag, role);
						}
					}
				}
			}	
		}

		// for(auto r : role) {
		// 	if(r == SLAVE) {
		// 		std::cout << "slave\n";
		// 	} else if(r == MASTER) {
		// 		std::cout << "master\n";
		// 	} else if(r == UNASSIGNED) {
		// 		std::cout << "unassigned\n";
		// 	}
		// }
	}

	template<class SpaceT, class FEBaseT>
	bool find_contacts(	SpaceT &space,
		const std::unique_ptr<FEBaseT> &master_fe, 
		const std::unique_ptr<FEBaseT> &slave_fe, 
		std::vector< std::shared_ptr<ContactAssembly> > &contacts,
		const libMesh::Real search_radius, const bool strict_gap_policy,
		const std::shared_ptr<moonolith::Predicate> &predicate)
	{
		using namespace libMesh;
		using namespace std;


		typedef Intersector::Scalar Scalar;
		static const Scalar tol = 1e-8;

		const MeshBase &mesh = space.mesh();
		const int dim = mesh.mesh_dimension();

		std::cout << "n_vol_elements:" << mesh.n_active_local_elem() << std::endl;

		// Chrono c;
		// c.start();
		std::vector<int> pairs;
		if(!boundary_hash_grid_detect_intersections(mesh, pairs, search_radius)) {
			std::cout << "NO INTERSECTIONS" << std::endl;
			return false;
		}

		// c.stop();
		std::cout << "proximity detection ";
		// c.describe(std::cout);


//		size_t n = mesh.n_active_local_elem();
		const int approx_order = space.order();

			 //init extra quadrature quantities
		slave_fe->get_xyz();


		Polyhedron polyhedron_1, polyhedron_2;
		DenseMatrix<Real> polygon_1, polygon_2, ref_polygon_1, ref_polygon_2, clipper_1;
		DenseMatrix<Real> side_polygon_1, side_polygon_2;
		DenseMatrix<Real> isect_polygon_1, isect_polygon_2;
//		Scalar isect_1[MAX_N_ISECT_POINTS * MAX_N_DIMS], isect_2[MAX_N_ISECT_POINTS * MAX_N_DIMS];

		DenseMatrix<Real> A(dim, dim), Ainv(dim, dim);
		DenseMatrix<Real> b(dim, 1), binv(dim, 1);

		shared_ptr<Transform> transform_1, transform_2;

		Intersector isector;

//		Scalar projection_1[MAX_N_ISECT_POINTS * MAX_N_DIMS];
//		Scalar projection_2[MAX_N_ISECT_POINTS * MAX_N_DIMS];

		std::vector<Point> points_1, points_2;


		Box box_1(dim), box_2(dim);

		QMortar q_1(dim), q_2(dim);
		QMortar qref_1(dim), qref_2(dim);


		size_t n_projections = 0;
		size_t n_candidates  = 0;

		Point n1, n2;

		bool intersected = false;
		Real local_element_matrices_sum = 0;

		contacts.reserve(pairs.size());

		std::shared_ptr<ContactAssembly> current_contact;

		for(auto it = pairs.begin(); it != pairs.end(); /*inside*/) {
		const auto index_1  = *it++;
		const auto index_2  = *it++;

		const Elem &el_1 = *mesh.elem(index_1);
		const Elem &el_2 = *mesh.elem(index_2);

		

			//FIXME This is a hack
		if(has_constrained_dofs(space, el_1) || 
			has_constrained_dofs(space, el_2)) {
			continue;
	}

	assert(is_valid_elem_type(el_1.type()));
	assert(is_valid_elem_type(el_2.type()));

	if(dim == 2) {
		make_polygon(el_1, polygon_1);
		make_polygon(el_2, polygon_2);

		transform_1 = make_shared<Transform2>(el_1);
		transform_2 = make_shared<Transform2>(el_2);

	} else if(dim == 3) {
		make_polyhedron(el_1, polyhedron_1);
		make_polyhedron(el_2, polyhedron_2);

		transform_1 = make_shared<Transform3>(el_1);
		transform_2 = make_shared<Transform3>(el_2);
	}

	for(uint side_1 = 0; side_1 < el_1.n_sides(); ++side_1) {
		if(el_1.neighbor_ptr(side_1) != nullptr) continue;
		auto side_ptr_1 = el_1.build_side_ptr(side_1);	

		compute_side_normal(dim, *side_ptr_1, n1);

		box_1.reset();
		enlarge_box_from_side(dim, *side_ptr_1, box_1, search_radius);

		if(dim == 2) {
			make_polygon(*side_ptr_1, side_polygon_1);
		} else if(dim == 3) {
			make_polygon_3(*side_ptr_1, side_polygon_1);
		} else {
			assert(false);
		}

		for(uint side_2 = 0; side_2 < el_2.n_sides(); ++side_2) {
			if(el_2.neighbor_ptr(side_2) != nullptr) continue;

			if(predicate) {
				if(!predicate->are_master_and_slave( 
					mesh.get_boundary_info().boundary_id(&el_1, side_1),
					mesh.get_boundary_info().boundary_id(&el_2, side_2)
					)) {
					continue;
			}
		}


		auto side_ptr_2 = el_2.build_side_ptr(side_2);	
		compute_side_normal(dim, *side_ptr_2, n2);

		const Real cos_angle = n1.contract(n2);

				//if the angle is more than 60 degrees ( cos(60/180*pi) == 0.5 ) or has same orientation skip
		if(cos_angle >= -0.5) {
			continue;
		}

		box_2.reset();
		enlarge_box_from_side(dim, *side_ptr_2, box_2, search_radius);

		if(!box_1.intersects(box_2, tol)) {
			continue;
		}



		++n_candidates;

		bool pair_intersected = false;

		if(dim == 2) {
			make_polygon(*side_ptr_2, side_polygon_2);

			if(!project_2D(side_polygon_1, side_polygon_2, isect_polygon_1, isect_polygon_2)) {
				continue;
			}

			const Scalar dx = side_polygon_2(0, 0) - side_polygon_2(1, 0);
			const Scalar dy = side_polygon_2(0, 1) - side_polygon_2(1, 1);

			const Scalar isect_dx = isect_polygon_2(0, 0) - isect_polygon_2(1, 0);
			const Scalar isect_dy = isect_polygon_2(0, 1) - isect_polygon_2(1, 1);

			const Scalar area   = std::sqrt(isect_dx*isect_dx + isect_dy*isect_dy);
			const Scalar area_slave = std::sqrt(dx*dx + dy*dy);
			const Scalar relative_area = area/area_slave;
			const Scalar weight = 1./area_slave;

			const int order = order_for_l2_integral(dim, el_1, approx_order, el_2, approx_order);

			make_composite_quadrature_on_surf_2D(isect_polygon_1, weight, order, q_1);
			make_composite_quadrature_on_surf_2D(isect_polygon_2, weight, order, q_2);

			pair_intersected = true;
			++n_projections;


			current_contact = std::make_shared<ContactAssembly>();
			current_contact->isect_area	   = area;
			current_contact->relative_area = relative_area;


		} else if(dim == 3) {
			make_polygon_3(*side_ptr_2, side_polygon_2);

			if(!project_3D(
				side_polygon_1, 
				side_polygon_2, 
				isect_polygon_1,
				isect_polygon_2))
			{
				continue;
			}

			const Scalar area_slave = isector.polygon_area_3(side_polygon_2.m(),  &side_polygon_2.get_values()[0]);	
			const Scalar area   	= isector.polygon_area_3(isect_polygon_2.m(), &isect_polygon_2.get_values()[0]);
			const Scalar relative_area 	= area/area_slave;
			const Scalar weight = 1./area_slave;

			const int order = order_for_l2_integral(dim, el_1, approx_order, el_2, approx_order);

			make_composite_quadrature_on_surf_3D(isect_polygon_1, weight, order, q_1);
			make_composite_quadrature_on_surf_3D(isect_polygon_2, weight, order, q_2);

			pair_intersected = true;
			++n_projections;


			current_contact = std::make_shared<ContactAssembly>();
			current_contact->isect_area	   = area;
			current_contact->relative_area = relative_area;

		} else {
			assert(false);
			return false;
		}

		if(pair_intersected) {
			transform_to_reference_surf(*transform_1, el_1.type(), q_1, qref_1);
			transform_to_reference_surf(*transform_2, el_2.type(), q_2, qref_2);

					// std::cout << "----------------------\n";
					// std::cout << "----------------------\n";
					// q_2.print_info();
					// std::cout << "----------------------\n";


			master_fe->attach_quadrature_rule(&qref_1);
					master_fe->reinit(&el_1);//, side_1);

					slave_fe->attach_quadrature_rule(&qref_2);
					slave_fe->reinit(&el_2);//, side_2);

					current_contact->parent_element_master = index_1;
					current_contact->side_number_master    = side_1;
					current_contact->id_master 			   = el_1.id();

					current_contact->parent_element_slave  = index_2;
					current_contact->side_number_slave 	   = side_2;
					current_contact->id_slave 			   = el_2.id();

					current_contact->coupling.zero();
					current_contact->gap.zero();
					current_contact->normals.zero();

				bool use_biorth_ = false; //ugly but works
					// bool use_biorth_ = true;
					libMesh::DenseMatrix<libMesh::Real> biorth_weights;
					
					if(use_biorth_) {

						std::unique_ptr<libMesh::FEVectorBase> biorth_elem = 
						libMesh::FEVectorBase::build(dim, space.dof_map().variable_type(space.var_num()));

						QMortar q_biorth(dim);

						const int order = order_for_l2_integral(dim, el_2, approx_order, el_2, approx_order);
						
						if(dim == 3) {
							make_composite_quadrature_on_surf_3D(side_polygon_2, 1., order, q_biorth);
						} else {
							make_composite_quadrature_on_surf_2D(side_polygon_2, 1., order, q_biorth);
						}

						QMortar q_biorth_ref(dim);
						transform_to_reference_surf(*transform_2, el_2.type(), q_biorth, q_biorth_ref);

						biorth_elem->attach_quadrature_rule(&q_biorth_ref);
						biorth_elem->reinit(&el_2);

						mortar_assemble_weights(*biorth_elem, biorth_weights);

						mortar_assemble_weighted_biorth(
							*master_fe, *slave_fe,
							biorth_weights,
							current_contact->coupling);


					} else {
						mortar_assemble(*master_fe, *slave_fe, current_contact->coupling);
					}

					const libMesh::Point pp = side_ptr_1->point(0);
					const Real plane_offset = n1.contract(pp);

					// std::cout << "-------------------------------------------------------------\n";
					// std::cout << "p_p = [" << pp(0) << ", " << pp(1) << ", " << pp(2) << "];\n";
					// std::cout << "n_p = [" << n1(0) << ", " << n1(1) << ", " << n1(2) << "];\n";
					// std::cout << "n_s = [" << n2(0) << ", " << n2(1) << ", " << n2(2) << "];\n";
					// std::cout << "-------------------------------------------------------------\n";
					// side_polygon_1.print();
					// std::cout << "\n=============================================================\n";
					// side_polygon_2.print();
					// std::cout << "\n=============================================================\n";

					// plot_polygon(dim, side_polygon_2.m(),  &side_polygon_2.get_values()[0],  "s/poly/"  + std::to_string(index_2));
					// plot_polygon(dim, isect_polygon_2.m(), &isect_polygon_2.get_values()[0], "s/isect/" + std::to_string(index_2));
					// plot_polygon(dim, side_polygon_1.m(),  &side_polygon_1.get_values()[0],  "m/poly/p" + std::to_string(index_1));
					// plot_polygon(dim, isect_polygon_1.m(), &isect_polygon_1.get_values()[0], "m/isect/" + std::to_string(index_1));
					// plot_quad_points(dim, q_2.get_points(), "slave/quadpoints");

					if(use_biorth_) {
						mortar_normal_and_gap_assemble_weighted_biorth(
							*slave_fe, 
							dim,
							n2,
							n1,
							plane_offset,
							biorth_weights,
							current_contact->normals, 
							current_contact->gap);

					} else {
						mortar_normal_and_gap_assemble(
							dim,
							*slave_fe, 
							n2,
							n1,
							plane_offset,
							current_contact->normals, 
							current_contact->gap);
					}

					const Scalar local_mat_sum = std::accumulate(current_contact->coupling.get_values().begin(), current_contact->coupling.get_values().end(), libMesh::Real(0.0));
					local_element_matrices_sum += local_mat_sum;
					intersected = true;

				// current_contact->describe();

                // std::cout<< "surface_assemble->isect_area = " << current_contact->isect_area <<std::endl;

                // std::cout<<" pow(surface_assemble->isect_area, dim/(dim-1.)) * dim = " << pow(current_contact->isect_area, dim/(dim-1.)) * dim  <<std::endl;

                // std::cout<<" partial sum = " << local_mat_sum <<std::endl;

                // std::cout<<" local_element_matrices_sum = " << local_element_matrices_sum <<std::endl;

                // assert(fabs(local_mat_sum - pow(current_contact->isect_area, dim/(dim-1.)) * dim) < 1e-8 || (!is_quad(el_2.type()) && !is_hex(el_2.type())));


					current_contact->finalize();


				// std::cout << "avg_gap: " << current_contact->avg_gap << std::endl;

					if(strict_gap_policy) {
						if(std::abs(current_contact->avg_gap) <= search_radius) {
							contacts.push_back(current_contact);	
							// plot_box(box_2, "intersected/box");
							// plot_polygon(dim, side_polygon_2.m(),  &side_polygon_2.get_values()[0],  "s/poly/"  + std::to_string(index_2));
						}
					} else {
						contacts.push_back(current_contact);
					}

					// abort();
				}
			}
		}
	}

	std::cout << "n_candidates:  " << n_candidates  << std::endl;
	std::cout << "n_projections: " << n_projections << std::endl; 
	std::cout << "local_element_matrices_sum: " << local_element_matrices_sum << std::endl;

	return intersected;
}

template<class SpaceT, class FEBaseT>
bool assemble_aux(
	SpaceT &space,
	const std::unique_ptr<FEBaseT> &master_fe, 
	const std::unique_ptr<FEBaseT> &slave_fe, 
	DSMatrixd &coupling, DVectord &gap, DVectord &normals, DSMatrixd &orthogonal_trafos, 
	std::vector<bool> &is_contact_node, const libMesh::Real search_radius, 
	const bool strict_gap_policy, 
	const std::shared_ptr<moonolith::Predicate> &predicate)
{
	using namespace libMesh;
	using namespace std;


	typedef Intersector::Scalar Scalar;
	static const Scalar tol = 1e-8;

	Intersector isector;

	const MeshBase &mesh = space.mesh();
	const int dim = mesh.mesh_dimension();

	std::vector< std::shared_ptr<ContactAssembly> > contacts;
	if(!find_contacts(space, master_fe, slave_fe, contacts, search_radius, strict_gap_policy, predicate)) return false;

	
	std::vector< std::vector<long> > adj_list;
	std::vector< std::vector<long> > dag;
	std::vector<long> ordering;
	
	build_adj_list(mesh, adj_list);
	build_dag(contacts, dag, ordering);


	std::vector<ushort> role;
	if(predicate) {
		role.resize(mesh.n_active_local_elem());
		std::fill(role.begin(), role.end(), SLAVE);

	} else {
		
		assign_master_and_slave_roles(dag, ordering, adj_list, role);
	}

	size_t n_dofs = space.dof_map().n_dofs();

	std::cout<<"*********************************************DOF = "<<n_dofs<<std::endl;

	gap = zeros(n_dofs);
	if(space.is_vector()) {
		normals = zeros(n_dofs);
	} else {
		normals = zeros(n_dofs * dim);
	}

	coupling = sparse(n_dofs, n_dofs, std::max(2, int(n_dofs * 0.1)));


	is_contact_node.resize(n_dofs); 
	std::fill(is_contact_node.begin(), is_contact_node.end(), false);

	{
		Write<DVectord>  w_gap(gap);
		Write<DSMatrixd> w_coupling(coupling);
		Write<DVectord>  w_normals(normals);

		std::vector<dof_id_type> dof_indices_slave;
		std::vector<dof_id_type> dof_indices_master;
		std::vector<dof_id_type> dof_indices_normal(dim, 0);

		for(dof_id_type dim_d = 0; dim_d < dim; ++dim_d) {
			dof_indices_normal[dim_d] = dim_d;
		}

		{
			//For performance reason add 0 entries to the diagonal when assemblying
			Range rr = row_range(coupling);
			for(SizeType i = rr.begin(); i < rr.end(); ++i) {
				coupling.add(i, i, 0.0);
			}
		}

		for(auto c_ptr : contacts) {
			//if the slave role is ok
			if(!c_ptr->is_valid || role[c_ptr->id_slave] != SLAVE) continue;
			auto ptr_slave  = mesh.elem(c_ptr->parent_element_slave);
			auto ptr_master = mesh.elem(c_ptr->parent_element_master);

			space.dof_map().dof_indices(ptr_slave,  dof_indices_slave,  space.var_num());
			space.dof_map().dof_indices(ptr_master, dof_indices_master, space.var_num());

			// c_ptr->describe(std::cout);
			// print_vector(dof_indices_slave.begin(), dof_indices_slave.end());
			// std::cout << "--------------------------\n" << std::endl;

			add_vector(c_ptr->gap, dof_indices_slave, gap);
			add_matrix(c_ptr->coupling, dof_indices_slave, dof_indices_master, coupling);

			DenseVector<Real> n_vec(c_ptr->normals.m()*c_ptr->normals.n());
			n_vec.get_values() = c_ptr->normals.get_values();


			// { //plot slave
			// 	auto side_ptr = ptr_slave->build_side_ptr( c_ptr->side_number_slave );	
			// 	DenseMatrix<Real>  side_polygon_2;
			// 	if(dim == 2) {
			// 		make_polygon(*side_ptr, side_polygon_2);
			// 	} else if(dim == 3) {
			// 		make_polygon_3(*side_ptr, side_polygon_2);
			// 	} else {
			// 		assert(false);
			// 	}

			// 	plot_polygon(dim, side_polygon_2.m(), &side_polygon_2.get_values()[0], "s/el/" + std::to_string(c_ptr->parent_element_slave));
			// }

			if(space.is_vector()) {
				std::vector<dof_id_type> side_dofs;
				side_dofs.reserve(ptr_slave->n_nodes());
				
				for(uint i = 0; i < ptr_slave->n_nodes(); ++i) {
					
					if(ptr_slave->is_node_on_side(i, c_ptr->side_number_slave)) {
						int sys_num = space.dof_map().sys_number();
						int var_num = space.var_num();

						for(unsigned int c = 0; c != ptr_slave->node_ref(i).n_comp(sys_num, var_num); ++c) {
							side_dofs.push_back(ptr_slave->node_ref(i).dof_number(sys_num, var_num,c));
						}
					} 
				}

				add_vector(n_vec, dof_indices_slave, normals);

				assert(side_dofs.size() == 4 || dim == 3);

				for(auto dof : side_dofs) {
					is_contact_node[dof] = true;
				}

			} else {
				auto ptr_side_slave  = ptr_slave->build_side_ptr(c_ptr->side_number_slave);

				std::vector<dof_id_type> dof_indices_slave_vec(dof_indices_slave.size() * dim, 0);
				
				for(uint i = 0; i < dof_indices_slave.size(); ++i) {
					for(uint d = 0; d < dim; ++d) {
						dof_indices_slave_vec[i*dim + d] = dof_indices_slave[i] * dim + d;
					}
				}

				add_vector(n_vec, dof_indices_slave_vec, normals);

				space.dof_map().dof_indices(ptr_side_slave.get(),  dof_indices_slave,  space.var_num());

				for(auto dof : dof_indices_slave) {
					is_contact_node[dof] = true;
				}
			}
		}
	}

	// write("coupling.m", coupling);
	// write("gap.m", gap);


	DVectord sum_c   = sum(coupling, 1);
	DVectord gap_2   = zeros(size(gap));


	// write("sum_c.m", sum_c);

	{
		Read<DVectord>  r_c(sum_c);
		Read<DVectord>  r_g(gap);
		Write<DVectord> w_g(gap_2);

		Range r = range(gap);
		for(uint i = r.begin(); i < r.end(); ++i) {
			const double val = sum_c.get(i);
			assert(std::abs(val) < 1e-8 ||  is_contact_node[i]);

			if(is_contact_node[i]) {
				gap_2.set(i, gap.get(i)/val);
			} 
		}
	}

	std::cout << "sum(coupling): " << double(sum(sum_c)) << std::endl;

	DSMatrixd coupling_2 = coupling;
	
	{
		Write<DSMatrixd> w_c(coupling_2);
		Read<DVectord>  r_c(sum_c);

		each_read(coupling, [&sum_c, &coupling_2, &is_contact_node](const SizeType i, const SizeType j, const Scalar value) {
			if(is_contact_node[i]) {
				coupling_2.set(i, j, value/sum_c.get(i));
			}
		});
	}

	gap = std::move(gap_2);
	coupling = std::move(coupling_2);


	DVectord mm_normals = normals;
	DVectord lenghts = zeros(n_dofs/dim);

	{
		Write<DVectord> w_l(lenghts);
		each_read(mm_normals, [&lenghts, dim](const int i, const double value){
			lenghts.add(i/dim, value*value);
		});
	}

	lenghts = sqrt(lenghts);

	{	
		Read<DVectord>  r_n(mm_normals);
		Read<DVectord>  r_l(lenghts);
		Write<DVectord> w_n(normals);

		Range r = range(normals);
		for(uint i = r.begin(); i < r.end(); ++i) {
			if(!is_contact_node[i]) {
				normals.set(i, 0.0);
			} else if(lenghts.get(i/dim)>1e-16) {
				normals.set(i, mm_normals.get(i)/lenghts.get(i/dim));
			}
		}
	}

	orthogonal_trafos = sparse(n_dofs, n_dofs, dim);

	{
		std::vector<Scalar> normal(dim, 0);
		std::vector<Scalar> H(dim * dim, 0);

		Read<DVectord>  r_n(normals);
		Write<DSMatrixd> w_o(orthogonal_trafos);

		Range r = range(normals);
		for(uint i = r.begin(); i < r.end(); i += dim) {
			bool use_identity = true;
			
			if(is_contact_node[i]) {
				
				for(uint d = 0; d < dim; ++d) {
					normal[d] = normals.get(i + d);
				}

				if(std::abs(normal[0] - 1.) > 1e-8) {
					use_identity = false;
					
					//-e1 basis vector
					normal[0] -= 1;

					Real len = 0.;
					for(uint d = 0; d < dim; ++d) {
						len += normal[d] * normal[d];
					}

					len = std::sqrt(len);

					assert(len > 0);

					for(uint d = 0; d < dim; ++d) {
						normal[d] /= len;
					}

					if(dim == 2) {
						isector.householder_reflection_2(&normal[0], &H[0]);
					} else {
						isector.householder_reflection_3(&normal[0], &H[0]);
					}


					for(uint di = 0; di < dim; ++di) {
						for(uint dj = 0; dj < dim; ++dj) {
							orthogonal_trafos.set(i + di, i + dj, H[di * dim + dj]);
						}
					}
				}
			} 

			if(use_identity)
			{
				for(uint di = 0; di < dim; ++di) {
					orthogonal_trafos.set(i + di, i + di, 1.);
				}
			}
		}
	}

	{
		Write<DVectord> w_g(gap);
		Range r = range(gap);

		const Scalar large_number = 1000;

		for(uint i = r.begin(); i < r.end(); ++i) {
			if(!is_contact_node[i] || i % dim > 0) {
				gap.set(i, large_number);
			}
		}
	}


	coupling += identity(size(coupling));
	return true;
}

bool MortarContactAssembler::assemble(DSMatrixd &coupling, DVectord &gap, DVectord &normals, DSMatrixd &orthogonal_trafos, std::vector<bool> &is_contact_node, const libMesh::Real search_radius, const std::shared_ptr<moonolith::Predicate> &predicate) {
	

	if(space_->is_vector()) {
		std::unique_ptr<libMesh::FEVectorBase> master_fe, slave_fe;
		master_fe = libMesh::FEVectorBase::build(space_->mesh().mesh_dimension(), 
			space_->dof_map().variable_type(space_->var_num()));

		slave_fe  = libMesh::FEVectorBase::build(space_->mesh().mesh_dimension(), 
			space_->dof_map().variable_type(space_->var_num()));
		return assemble_aux(*space_, master_fe, slave_fe, coupling, gap, normals, orthogonal_trafos, is_contact_node, search_radius, strict_gap_policy, predicate);

	} else {
		std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
		master_fe = libMesh::FEBase::build(space_->mesh().mesh_dimension(), 
			space_->dof_map().variable_type(space_->var_num()));

		slave_fe  = libMesh::FEBase::build(space_->mesh().mesh_dimension(), 
			space_->dof_map().variable_type(space_->var_num()));
		return assemble_aux(*space_, master_fe, slave_fe, coupling, gap, normals, orthogonal_trafos, is_contact_node, search_radius, strict_gap_policy, predicate);
	}
}
}

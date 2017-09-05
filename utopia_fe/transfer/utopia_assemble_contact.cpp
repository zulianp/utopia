#include "utopia_assemble_contact.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_LibMeshBackend.hpp"

#include "MortarAssemble.hpp"
#include "MortarAssembler.hpp"

#include "utopia_STree.hpp"
#include "utopia_ElementDofMap.hpp"
#include "utopia_FESpaceAdapter.hpp"

// #include "MapSparseMatrix.hpp"

#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/serial_mesh.h"

#include "moonolith_profiler.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_sparse_matrix.hpp"

#include "utopia_Socket.hpp"

#include <cmath>
#include <queue>

namespace utopia {

	using namespace libMesh;
	
	template<typename T>
	static void normalize(std::vector<T> &vec)
	{
		T len = 0.;
		for(uint d = 0; d < vec.size(); ++d) {
			len += vec[d] * vec[d];
		}
		
		len = std::sqrt(len);
		
		assert(len > 0);
		
		for(uint d = 0; d < vec.size(); ++d) {
			vec[d] /= len;
		}
	}
	
	template<typename T>
	inline void Print(const std::vector<T> &v, std::ostream &os)
	{
		for(auto i : v) {
			os << i << " ";
		}
		
		os << "\n";
	}
	
	template<int Dimensions>
	inline static void nodes_are_boundary_hack(const libMesh::Elem &el,
											   const int n_dofs,
											   const int approx_order,
											   const int side,
											   std::vector<bool> &result)
	{
	
		result.resize(n_dofs, false);

		std::vector<unsigned int> dofs;
		libMesh::FE<Dimensions, LAGRANGE>::dofs_on_side(&el, libMesh::Order(approx_order), 0, dofs);
		
		for(auto d : dofs) {
			assert(d < result.size());
			result[d] = true;
		}
	}


	inline static void nodes_are_boundary_hack(const libMesh::FEBase &fe,
											   std::vector<bool> &result)
	{
		auto &f = fe.get_phi();
		result.resize(f.size(), false);

		for(std::size_t i = 0; i < f.size(); ++i) {
			for(std::size_t qp = 0; qp < f[i].size(); ++qp) {
				if(f[i][qp] > 1e-16) {
					result[i] = true;
				}
			}
		}
	}

	inline static bool check_boundary_and_fun_consistent(const libMesh::FEBase &fe, const std::vector<bool> &node_is_boundary)
	{
		auto &f = fe.get_phi();
		for(std::size_t i = 0; i < node_is_boundary.size(); ++i) {
			for(std::size_t qp = 0; qp < f[i].size(); ++qp) {
				if(f[i][qp] > 1e-8) {
					assert(node_is_boundary[i]);
					if(!node_is_boundary[i]) return false;
				}
			}
		}

		return true;
	}

	inline static bool check_positive_funs(libMesh::FEBase &slave_fe)
	{
		auto &f = slave_fe.get_phi();

		for(auto &f_i : f) {
			for(auto f_iq : f_i) {
				if(f_iq < -1e-8) {
					assert(f_iq >= -1e8);
					return false;
				}
			}
		}

		return true;
	}
	
	inline static bool check_node_is_boundary(const ElemType &type,
											  const std::vector<bool> &is_boundary)
	{
		
		int n_bound = 0;
		for(std::size_t i = 0; i < is_boundary.size(); ++i) {
			n_bound += is_boundary[i];
		}
		
		switch(type) {
			case TRI3:
			case TRISHELL3:
			case QUAD4:
			case QUADSHELL4:
			{
				assert(n_bound == 2);
				return n_bound == 2;
			}
				
			case TET4:
			{
				assert(n_bound == 3);
				return n_bound == 3;
			}
				
			default:
			{
				std::cerr << "Missing implementation for ElemType " << type << std::endl;
				assert(false && "implement me!");
				break;
			}
		}
		
		return false;
	}
	
	//	static bool check_lumped_is_positive(const libMesh::DenseMatrix<libMesh::Real> &mat)
	//	{
	//		std::vector<libMesh::Real> lumped(mat.m(), 0.);
	//
	//		for(int i = 0; i < mat.m(); ++i) {
	//			for(int j = 0; j < mat.n(); ++j) {
	//				lumped[i] += mat(i, j);
	//			}
	//		}
	//
	//
	//		for(auto v : lumped) {
	//			assert(v >= 0.);
	//			if(v < 0.) return false;
	//		}
	//
	//		return true;
	//	}
	//
	inline static void assemble_trace_biorth_weights_from_space(const ElemType &type,
																const std::vector<bool> &is_boundary,
																libMesh::DenseMatrix<libMesh::Real> &weights)
	{
		switch(type) {
			case TRI3:
			case TRISHELL3:
			{
				weights.resize(3, 3);
				weights.zero();
				
#ifndef NDEBUG
				int n_bound = 0;
				for(std::size_t i = 0; i < is_boundary.size(); ++i) {
					n_bound += is_boundary[i];
				}
				
				assert(n_bound == 2);
#endif //NDEBUG
				
				for(std::size_t i = 0; i < is_boundary.size(); ++i) {
					if(!is_boundary[i]) { continue; }
					
					weights(i, i) = 2;
					
					for(std::size_t j = 0; j < is_boundary.size(); ++j) {
						if(is_boundary[j] && i != j) {
							weights(i, j) = -1;
						}
					}
				}
				
				break;
			}
				
			case QUAD4:
			case QUADSHELL4:
			{
				
				weights.resize(4, 4);
				weights.zero();
#ifndef NDEBUG
				int n_bound = 0;
				for(std::size_t i = 0; i < is_boundary.size(); ++i) {
					n_bound += is_boundary[i];
				}
				
				assert(n_bound == 2);
#endif //NDEBUG
				
				for(std::size_t i = 0; i < is_boundary.size(); ++i) {
					if(!is_boundary[i]) { continue; }
					
					weights(i, i) = 2;
					
					for(std::size_t j = 0; j < is_boundary.size(); ++j) {
						if(is_boundary[j] && i != j) {
							weights(i, j) = -1;
						}
					}
				}
				
				break;
			}
				
			case TET4:
			{
				weights.resize(4, 4);
				weights.zero();
				
#ifndef NDEBUG
				int n_bound = 0;
				for(std::size_t i = 0; i < is_boundary.size(); ++i) {
					n_bound += is_boundary[i];
				}
				
				assert(n_bound == 3);
#endif //NDEBUG
				
				for(std::size_t i = 0; i < is_boundary.size(); ++i) {
					if(!is_boundary[i]) { continue; }
					
					weights(i, i) = 3;
					
					for(std::size_t j = 0; j < is_boundary.size(); ++j) {
						if(is_boundary[j] && i != j) {
							weights(i, j) = -1;
						}
					}
				}
				
				break;
			}
				
			default: {
				std::cerr << "Missing implementation for ElemType " << type << std::endl;
				assert(false && "implement me!");
				break;
			}
		}
	}
	
	template<class Iterator>
	static void write_space(
							const Iterator &begin,
							const Iterator &end,
							MeshBase &space,
							const std::vector<ElementDofMap> &dof_map,
							const std::vector<long> &variable_number,
							const std::vector<long> &variable_order,
							const std::vector<ElementDofMap> &subdomain_id,
							const std::vector<ElementDofMap> &side_set_id,
							const std::vector<ElementDofMap> &face_set_id_global,
							moonolith::OutputStream &os)
	{
		const int dim 		  = space.mesh_dimension();
		const long n_elements = std::distance(begin, end);
		
		std::set<long> nodeIds;
		std::map<long, long> mapping;
		std::vector<dof_id_type> dof_array;
		
		for(Iterator it = begin; it != end; ++it) {
			
			const Elem *elem = space.elem(*it);
			
			for(dof_id_type j = 0; j != elem->n_nodes(); ++j) {
				
				nodeIds.insert(elem->node(j));
				
				
			}
		}
		
		long n_nodes = nodeIds.size();
		
		// Estimate for allocation
		os.request_space( (n_elements * 8 + n_nodes * dim) * (sizeof(double) + sizeof(long)) );
		
		//WRITE 1
		os << dim;
		
		int index = 0;
		for (auto nodeId : nodeIds) {
			mapping[nodeId] = index++;
		}
		
		//WRITE 2
		os << n_nodes;
		
		//WRITE 6
		os << n_elements;
		
		for(auto node_id : nodeIds){
			
			const Point &p = space.node(node_id);
			
			for(int i = 0; i < dim; ++i) {
				
				//WRITE 3
				os << p(i);
				
			}
		}
		
		std::vector<dof_id_type> indices_vector;
		
		
		
		for(Iterator it = begin; it != end; ++it) {
			
			//			const int k = *it;
			
			const Elem *elem = space.elem(*it);
			
			const int e_n_nodes = elem->n_nodes();
			
			const int type = elem->type();
			
			//WRITE 7
			os << type << e_n_nodes;
			
			
			
			for (int i = 0; i != e_n_nodes; ++i) {
				
				auto it = mapping.find(elem->node(i));
				
				assert(it != mapping.end());
				
				int index = it->second;
				
				//WRITE 8
				os << index;
				
			}
			
			
			
			//WRITE 9
			assert(!dof_map.at(elem->id()).empty());
			
			os << dof_map.at(elem->id());
			
			//			bool  size=true;
			
			int volume_tag;
			
			volume_tag=subdomain_id[elem->id()].global.at(0);
			
			os << volume_tag;
			
			int side_set_tag;
			
			//			int face_id;
			
			//			bool check_side_id_one = true;
			side_set_tag=side_set_id[elem->id()].global.at(0);
			
			os << side_set_tag;
			os << face_set_id_global.at(elem->id());
		}
		//
		//
		
		//WRITE 11
		os << variable_number.at(0);
		
		//WRITE 12
		os << variable_order.at(0);
		
	}
	
	template<class Iterator>
	static void write_element_selection(
										const Iterator &begin,
										const Iterator &end,
										const FESpaceAdapter &fespace,
										moonolith::OutputStream &os)
	{
		write_space(begin, end, *fespace.mesh(), fespace.dof_map(), fespace.variable_number(), fespace.variable_number(), fespace.subdomain_id(), fespace.side_set_id(), fespace.face_set_id_global(), os);
	}
	
	
	static void read_space(moonolith::InputStream &is, std::shared_ptr<MeshBase> & space,
						   std::vector<ElementDofMap> &dof_map,
						   std::vector<long> &variable_number,
						   std::vector<long> &variable_order,
						   std::vector<ElementDofMap> &subdomain_id,
						   std::vector<ElementDofMap> &side_set_id,
						   std::vector<ElementDofMap> &face_set_id_global,
						   const libMesh::Parallel::Communicator &comm)
	{
		using namespace std;
		
		//READ 1
		int dim;
		is >> dim;
		
		//READ 2
		long n_nodes;
		is >> n_nodes;
		
		//READ 6
		long n_elements;
		is >> n_elements;
		
		auto mesh_ptr = std::make_shared<SerialMesh>(comm, dim);
		
		mesh_ptr->reserve_nodes(n_nodes);
		
		for (long iii = 0; iii != n_nodes; ++iii) {
			
			Point p;
			
			for(int j = 0; j < dim; ++j) {
				//READ 3
				is >> p(j);
			}
			
			mesh_ptr->add_point(p);
		}
		
		
		
		dof_map.resize(n_elements);
		
		subdomain_id.resize(n_elements);
		
		side_set_id.resize(n_elements);
		
		face_set_id_global.resize(n_elements);
		
		face_set_id_global.resize(n_elements);
		
		
		for(long i = 0; i !=n_elements; ++i) {
			
			//READ 7
			
			int type, e_n_nodes;
			
			is >> type >> e_n_nodes;
			
			//std::cout<<"e_n_nodes_read = "<<e_n_nodes<<std::endl;
			
			auto elem =  Elem::build(ElemType(type)).release();
			
			//std::cout<<"n_side_read ="<< elem->n_sides()<<std::endl;
			
			
			int index;
			
			for (int ii = 0; ii != e_n_nodes; ++ii) {
				
				//READ 8
				is >> index;
				//std::cout<<"index = "<<index<<std::endl;
				elem->set_node(ii) = & mesh_ptr->node(index);
				
			}
			
			
			//READ 9
			is >> dof_map.at(i);
			//std::cout<< "dof_map_read = "<<dof_map[i].global.at(0)<<std::endl;
			
			int volume_tag, side_set_tag;
			//			int face_id;
			
			//			bool on_boundary=false;
			//std::cout<<"read n_elements = "<<n_elements<<std::endl;
			
			
			is >> volume_tag;
			
			//std::cout<<" read volume role = "<< volume_tag <<std::endl;
			
			subdomain_id[i].global.insert(subdomain_id[i].global.end(),volume_tag);
			
			is >> side_set_tag;
			
			is >> face_set_id_global.at(i);
			
			//std::cout <<"read value"<< face_set_id_global[i].global.at(0)<<std::endl;
			
			side_set_id[i].global.insert(side_set_id[i].global.end(),side_set_tag);
			
			mesh_ptr->add_elem(elem);
			
			libmesh_assert(elem);
			
		}
		
		//READ 11
		variable_number.resize(1);
		is >> variable_number.at(0);
		
		//READ 12
		variable_order.resize(1);
		is >> variable_order.at(0);
		
		
		//!!!! dummy parameters
		space = mesh_ptr;
		
	}
	
	static void read_spaces(moonolith::InputStream &is, FESpaceAdapter &utopiamesh, const libMesh::Parallel::Communicator &comm_mesh)
	{
		read_space(is, utopiamesh.mesh(), utopiamesh.dof_map(), utopiamesh.variable_number(), utopiamesh.variable_order(), utopiamesh.subdomain_id(), utopiamesh.side_set_id(), utopiamesh.face_set_id_global(), comm_mesh);
	}
	
	
	template<int Dimensions, class Fun>
	static bool SurfaceAssemble(moonolith::Communicator &comm,
								const std::shared_ptr<MeshBase> &master_slave,
								const std::shared_ptr<DofMap> &dof_map,
								const unsigned int var_num,
								Fun process_fun,
								const moonolith::SearchSettings &settings,
								const libMesh::Real search_radius,
								const std::vector< std::pair<int, int> >  &tags,
								const bool use_biorth)
	{
		
		std::shared_ptr<FESpaceAdapter> local_fun_spaces_new = std::make_shared<FESpaceAdapter>(master_slave, dof_map, var_num, tags);
		auto predicate = std::make_shared<moonolith::MasterAndSlave>();
		
		for(auto t : tags)
			predicate->add(t.first, t.second);
		
		using namespace moonolith;

		
		typedef STree<Dimensions> NTreeT;
		typedef typename NTreeT::DataContainer DataContainer;
		typedef typename NTreeT::DataType SurfaceAdapter;
		
		long maxNElements = settings.max_elements;
		long maxDepth = settings.max_depth;
	
		
		const auto &mesh = master_slave;
		
		const int n_elements = mesh->n_elem();
		
		const Parallel::Communicator &libmesh_comm_mesh = master_slave->comm();
		
		//		const int dim_master = master_slave->mesh_dimension();
		//		const int dim_slave = master_slave->mesh_dimension();
		
		MeshBase::const_element_iterator e_it = mesh->active_elements_begin();
		const MeshBase::const_element_iterator e_end = mesh->active_elements_end();
		std::vector<int> block_id;
		std::vector<int> block_id_def;
		
		MOONOLITH_EVENT_BEGIN("create_adapters");
		////////////////////////////////////////////////////////////////////////////////////////////////////
		std::shared_ptr<NTreeT> tree = NTreeT::New(predicate, maxNElements, maxDepth);
		tree->reserve(n_elements);
				
		std::shared_ptr<FESpaceAdapter> local_spaces = std::make_shared<FESpaceAdapter>(master_slave, dof_map, var_num, tags);
				
		for (auto it = master_slave->active_local_elements_begin();
			 it != master_slave->active_local_elements_end(); ++it) {
			
			auto elem = *it;
			
			if(!elem->on_boundary()) {
				continue;
			}
			
			bool check_size=false;
			
			for(uint side_elem = 0; side_elem < elem->n_sides(); ++side_elem){
				if ((predicate->select(master_slave->get_boundary_info().boundary_id(elem, side_elem))) && check_size==false){
					SurfaceAdapter a(*master_slave, elem->id(), elem->id(), master_slave->get_boundary_info().boundary_id(elem, side_elem), search_radius);
					assert(!local_spaces->dof_map()[elem->id()].empty());
					a.set_dof_map(&local_spaces->dof_map()[elem->id()].global);
					a.set_face_id(&local_spaces->face_set_id_global()[elem->id()].global);
					tree->insert(a);
					check_size=true;
				}
			}
		}
		
		tree->root()->bound().static_bound().enlarge(1e-8);
	
		
		////////////////////////////////////////////////////////////////////////////////////////////////////
		MOONOLITH_EVENT_END("create_adapters");
		
		//Just to have an indexed-storage
		std::map<long, std::shared_ptr<FESpaceAdapter> > utopiamesh;
		std::map<long, std::vector<std::shared_ptr<FESpaceAdapter> > > migrated_meshes;
		
		auto read = [&utopiamesh, &migrated_meshes, block_id, comm, &libmesh_comm_mesh, search_radius]
		(
		 const long ownerrank,
		 const long senderrank,
		 bool is_forwarding, DataContainer &data,
		 InputStream &in
		 ) {
			
			CHECK_STREAM_READ_BEGIN("vol_proj", in);
			
			auto proc_space = std::make_shared<FESpaceAdapter>(comm);
			
			read_spaces(in, *proc_space, libmesh_comm_mesh);
			
			if (!is_forwarding) {
				assert(!utopiamesh[ownerrank]);
				utopiamesh[ownerrank] = proc_space;
			} else {
				migrated_meshes[ownerrank].push_back(proc_space);
			}
			
			data.reserve(data.size() + proc_space->n_elements());
			
			auto s = proc_space->mesh();
			
			for (int i = 0; i<s->n_elem(); ++i) {
				int tag =proc_space->side_set_id()[i].global.at(0);
				data.push_back(SurfaceAdapter(*s, i, i,tag,search_radius));
				assert(!proc_space->dof_map()[i].empty());
				assert(!proc_space->side_set_id()[i].empty());
				data.back().set_dof_map(&proc_space->dof_map()[i].global);
				data.back().set_face_id(&proc_space->face_set_id_global()[i].global);
			}
			
			CHECK_STREAM_READ_END("vol_proj", in);
		};
		
		
		
		auto write = [&local_spaces, &utopiamesh, &comm]
		(
		 const long ownerrank, const long recvrank,
		 const std::vector<long>::const_iterator &begin,
		 const std::vector<long>::const_iterator &end,
		 const DataContainer &data,
		 OutputStream &out) {
			
			CHECK_STREAM_WRITE_BEGIN("vol_proj", out);
			
			if (ownerrank == comm.rank()) {
				write_element_selection(begin, end, *local_spaces, out);
				
				
			} else {
				auto it = utopiamesh.find(ownerrank);
				assert(it != utopiamesh.end());
				std::shared_ptr<FESpaceAdapter> spaceptr = it->second;
				assert(std::distance(begin, end) > 0);
				write_element_selection(begin, end, *spaceptr, out);
				
			}
			
			CHECK_STREAM_WRITE_END("vol_proj", out);
		};
		
		long n_false_positives = 0, n_projections = 0;

		auto fun = [&n_false_positives, &n_projections, &process_fun](
																	  SurfaceAdapter &master, SurfaceAdapter &slave) -> bool {
			bool ok = process_fun(master, slave);
			
			if(ok) {
				n_projections++;
				return true;
			} else {
				n_false_positives++;
				return false;
			}
			return true;
			
		};
				
		moonolith::search_and_compute(comm, tree, predicate, read, write, fun, settings);
		
		long n_total_candidates = n_projections + n_false_positives;
		
		long n_collection[3] = {n_projections, n_total_candidates, n_false_positives};
		comm.all_reduce(n_collection, 3, moonolith::MPISum());
		
		if (comm.is_root()) {
			std::cout << "n_intersections: " << n_collection[0]
			<< ", n_total_candidates: " 	 << n_collection[1]
			<< ", n_false_positives: " 	     << n_collection[2] << std::endl;
		}
		
		return true;
	}
	
	template<int Dimensions>
	bool SurfaceAssemble(
						 moonolith::Communicator &comm,
						 const std::shared_ptr<MeshBase> &master_slave,
						 const std::shared_ptr<DofMap> &dof_map,
						 const unsigned int var_num,
						 DSMatrixd &B,
						 DSMatrixd &orthogonal_trafos,
						 DVectord &gap,
						 DSMatrixd &normals,
						 DVectord &is_contact_node,
						 const moonolith::SearchSettings &settings,
						 const libMesh::Real search_radius,
						 const std::vector< std::pair<int, int> > &tags,
						 const bool use_biorth,
						 const bool use_volume_differential)
	{
		std::shared_ptr<FESpaceAdapter> local_fun_spaces_new = std::make_shared<FESpaceAdapter>(master_slave, dof_map, var_num, tags);
		
		libMesh::DenseMatrix<libMesh::Real> points_master;
		libMesh::DenseMatrix<libMesh::Real> points_slave;
		libMesh::DenseMatrix<libMesh::Real> intersection2;
		Polyhedron poly_master, poly_slave;
		Intersector isector;
		
		auto predicate = std::make_shared<moonolith::MasterAndSlave>();
		
		for(auto t : tags){
			predicate->add(t.first, t.second);
			std::cout << "[Status] Added master-slave pair = "<< t.first << ", " << t.second << std::endl;
		}
		
		static const double tol = 1e-8;
		
		std::vector<libMesh::dof_id_type> dofs_master, dofs_slave;
		
		libMesh::DenseMatrix<libMesh::Real> elemmat;//, massmat;
		libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;
		
		DenseMatrix<Real> side_polygon_master, side_polygon_slave;
		DenseMatrix<Real> isect_polygon_master, isect_polygon_slave;
		
		std::shared_ptr<Transform> trafo_master;
		std::shared_ptr<Transform> trafo_slave;
		
		Point n_master, n_slave;
		
		const int dim = master_slave->mesh_dimension();
		
		libMesh::Real local_element_matrices_sum = 0.0;
		
		//all face dof buffers
		moonolith::SparseMatrix<double> B_buffer(comm);
		moonolith::SparseMatrix<double> P_buffer(comm);
		moonolith::SparseMatrix<double> Q_buffer(comm);
		moonolith::SparseMatrix<double> rel_area_buffer(comm);
		moonolith::SparseMatrix<double> gap_buffer(comm);
		moonolith::SparseMatrix<double> normal_buffer(comm);
		
		DenseMatrix<Real> biorth_weights;
		
		auto fun = [&](const SElementAdapter<Dimensions> &master,
					   const SElementAdapter<Dimensions> &slave) -> bool {
			
			using namespace moonolith;

			const auto &master_mesh = master.space();
			const auto &slave_mesh  = slave.space();
			
			libMesh::DenseMatrix<libMesh::Real> elemmat;
			
			const int index_master  = master.element();
			const int index_slave = slave.element();
			
			auto &el_master  = *master_mesh.elem(index_master);
			auto &el_slave = *slave_mesh.elem(index_slave);
			
			const int dim_master = master_mesh.mesh_dimension();
			const int dim_slave = slave_mesh.mesh_dimension();
			
			Box box_master(dim_master), box_slave(dim_slave);
			
			QMortar ir_ref_master(dim_master);
			QMortar ir_master(dim_master);
			QMortar ir_slave(dim_slave);
			QMortar ir_ref_slave(dim_slave);
			
			//only works because there are not mixed elements
			const int approx_order = local_fun_spaces_new->variable_order()[0];
			
			std::shared_ptr<Contact> surface_assemble;
			
			const auto &side_id_master = master.dof_map_face();
			const auto &side_id_slave  = slave.dof_map_face();
			
			std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
			
			master_fe = libMesh::FEBase::build(master_mesh.mesh_dimension(), libMesh::Order(approx_order));
			slave_fe  = libMesh::FEBase::build(slave_mesh.mesh_dimension(),  libMesh::Order(approx_order));

		
			
			master_fe->get_phi();
//			master_fe->get_JxW(); //not necessary
			
			slave_fe->get_xyz();
			slave_fe->get_phi();
			slave_fe->get_JxW();


			//begin: hack stuff
			std::unique_ptr<libMesh::FEBase>  master_fe_hack, slave_fe_hack;
			master_fe_hack = libMesh::FEBase::build(master_mesh.mesh_dimension(), libMesh::Order(approx_order));
			slave_fe_hack  = libMesh::FEBase::build(slave_mesh.mesh_dimension(),  libMesh::Order(approx_order));
			
			libMesh::QGauss ir_hack(dim-1, libMesh::Order(1)); 
			
			if(dim == 2) {
				ir_hack.init(libMesh::EDGE2);
			} else {
				ir_hack.init(libMesh::TRI6);
			}

			master_fe_hack->get_phi();
			slave_fe_hack->get_phi(); 

			master_fe_hack->attach_quadrature_rule(&ir_hack);
			slave_fe_hack->attach_quadrature_rule(&ir_hack);

			//end: hack stuff

			typedef Intersector::Scalar Scalar;
			
			if(dim_slave == 2)  {
				make_polygon(el_master, points_master);
				make_polygon(el_slave,  points_slave);
				
				if(use_volume_differential) {
					trafo_master = std::make_shared<AffineTransform2>(el_master);
					trafo_slave  = std::make_shared<AffineTransform2>(el_slave);
				}
				
			} else if(dim_slave == 3) {
				make_polyhedron(el_master, poly_master);
				make_polyhedron(el_slave,  poly_slave);
				
				if(use_volume_differential) {
					trafo_master = std::make_shared<AffineTransform3>(el_master);
					trafo_slave  = std::make_shared<AffineTransform3>(el_slave);
				}
			}
			
			bool intersected = false;
			
			for(uint side_index_master = 0; side_index_master < el_master.n_sides(); ++side_index_master) {
				if(side_id_master[side_index_master] < 0) continue;
				
				if(el_master.neighbor_ptr(side_index_master) != nullptr) {
					std::cerr << "[Warning] it should never happen" << std::endl;
					continue;
				}
				
				auto side_master = el_master.build_side_ptr(side_index_master);
				
				compute_side_normal(dim_master, *side_master, n_master);
				
				if(fix_normal_orientation(el_master, side_index_master, n_master)) {
					std::cerr << "[Warning] fixed normal orientation of master face" << std::endl;
				}
				
				box_master.reset();
				enlarge_box_from_side(dim_master, *side_master, box_master, search_radius);
				
				if(dim_master == 2) {
					make_polygon(*side_master, side_polygon_master);
				} else if(dim_master == 3) {
					make_polygon_3(*side_master, side_polygon_master);
				} else {
					assert(false);
				}
				
				for(uint side_index_slave = 0; side_index_slave < el_slave.n_sides(); ++side_index_slave) {
					if(side_id_slave[side_index_slave] < 0) continue;
					
					if(el_slave.neighbor_ptr(side_index_slave) != nullptr) {
						std::cerr << "[Warning] it should never happen" << std::endl;
						continue;
					}
					//FIXME at some point we migth want to discrimate between faces of the same element
					// if (!predicate->tagsAreRelated(tag_1, tag_2)) continue;
					
					auto side_slave = el_slave.build_side_ptr(side_index_slave);
					
					compute_side_normal(dim_slave, *side_slave, n_slave);
					
					if(fix_normal_orientation(el_slave, side_index_slave, n_slave)) {
						std::cerr << "[Warning] fixed normal orientation of slave face" << std::endl;
					}
					
					const Real cos_angle = n_master.contract(n_slave);
					
					//if the angle is more than 60 degrees ( cos(60/180*pi) == 0.5 ) or has same orientation skip
					if(cos_angle >= -0.5) {
						continue;
					}
					
					box_slave.reset();
					enlarge_box_from_side(dim_slave, *side_slave, box_slave, search_radius);
					
					if(!box_master.intersects(box_slave, tol)) {
						continue;
					}
					
					bool pair_intersected = false;
					if(dim_slave == 2){
						make_polygon(*side_slave, side_polygon_slave);
						
						//plot_lines(2, 2, &side_polygon_master.get_values()[0], "in_master/" + std::to_string(master_facq) + "_" + std::to_string(cos_angle));
						//plot_lines(2, 2, &side_polygon_slave.get_values()[0], "in_slave/" + std::to_string(side_id_slave[0]) + "_" + std::to_string(cos_angle));
						
						if(!project_2D(side_polygon_master, side_polygon_slave, isect_polygon_master, isect_polygon_slave)){
							continue;
						}
						
						
						const Scalar dx = side_polygon_slave(0, 0) - side_polygon_slave(1, 0);
						const Scalar dy = side_polygon_slave(0, 1) - side_polygon_slave(1, 1);
						
						const Scalar isect_dx = isect_polygon_slave(0, 0) - isect_polygon_slave(1, 0);
						const Scalar isect_dy = isect_polygon_slave(0, 1) - isect_polygon_slave(1, 1);
						
						const Scalar area = std::sqrt(isect_dx*isect_dx + isect_dy*isect_dy);
						const Scalar area_slave = std::sqrt(dx*dx + dy*dy);
						const Scalar relative_area = area/area_slave;
						const Scalar weight = 1./area_slave;
						
						if(weight < 1e-15) continue;
						
						const int order = order_for_l2_integral(dim_master, el_master, approx_order, el_slave, approx_order);
						
						make_composite_quadrature_on_surf_2D(isect_polygon_master, weight, order, ir_master);
						make_composite_quadrature_on_surf_2D(isect_polygon_slave, weight, order, ir_slave);
						
						//						//Maybe remove
						//						ir_master.get_weights() = ir_slave.get_weights();
						
						pair_intersected = true;
						
						surface_assemble = std::make_shared<Contact>();
						surface_assemble->isect_area	   = area;
						surface_assemble->relative_area    = relative_area;
						
						// plot_polygon(2, 2, &side_polygon_master.get_values()[0], "master");
						// plot_polygon(2, 2, &side_polygon_slave.get_values()[0], "slave");
						
						
					} else if(dim_slave == 3) {
						make_polygon_3(*side_slave, side_polygon_slave);
						
						if(!project_3D(
									   side_polygon_master,
									   side_polygon_slave,
									   isect_polygon_master,
									   isect_polygon_slave)) {
							continue;
						}
						
						const Scalar area_slave = isector.polygon_area_3(side_polygon_slave.m(),  &side_polygon_slave.get_values()[0]);
						const Scalar area   	= isector.polygon_area_3(isect_polygon_slave.m(), &isect_polygon_slave.get_values()[0]);
						const Scalar relative_area = area/area_slave;
						const Scalar weight = 1./area_slave;
						
						assert(area_slave > 0);
						assert(area > 0);
						assert(weight > 0);
						
						const int order = order_for_l2_integral(dim_master, el_master, approx_order, el_slave, approx_order);
						
						make_composite_quadrature_on_surf_3D(isect_polygon_master, weight, order, ir_master);
						make_composite_quadrature_on_surf_3D(isect_polygon_slave,  weight, order, ir_slave);
						
						//Maybe remove
						ir_master.get_weights() = ir_slave.get_weights();
						
						pair_intersected = true;
						
						surface_assemble = std::make_shared<Contact>();
						surface_assemble->isect_area	= area;
						surface_assemble->relative_area = relative_area;
						
					} else {
						assert(false);
						return false;
					}
					
					// const bool enable_vis = side_id_slave.at(side_index_slave) == 434;
					// const bool enable_vis = true; //visdbg
					
					if(pair_intersected) {
						// if(enable_vis) { //visdbg
						// 	plot_polygon(dim, isect_polygon_master.get_values().size()/dim, &isect_polygon_master.get_values()[0], "master_" + std::to_string(side_id_slave[side_index_master])); //visdbg
						// 	plot_polygon(dim, isect_polygon_slave.get_values().size()/dim,  &isect_polygon_slave.get_values()[0],  "slave_" + std::to_string(side_id_slave[side_index_slave]));  //visdbg
						// } //visdbg
						
						// std::cout << "isect: " << master.handle() << " -> " << slave.handle() << std::endl;
						
						if(!use_volume_differential) {
							switch(dim) {
								case 2: {
									trafo_master = std::make_shared<SideAffineTransform2>(el_master, side_index_master);
									trafo_slave  = std::make_shared<SideAffineTransform2>(el_slave, side_index_slave);
									break;
								}
								case 3: {
									trafo_master = std::make_shared<SideAffineTransform3>(el_master, side_index_master);
									trafo_slave  = std::make_shared<SideAffineTransform3>(el_slave, side_index_slave);
									break;
								}
								default: {
									assert(false);
									break;
								}
							}
						}
						
						//////////////////////////////////ASSEMBLY ////////////////////////////////////////
						//////////////////////////////////////////////////////////////////////////////////////
						transform_to_reference_surf(*trafo_master, el_master.type(),  ir_master, ir_ref_master);
						transform_to_reference_surf(*trafo_slave,  el_slave.type(),   ir_slave,  ir_ref_slave);
						
						//master fe init
						master_fe->attach_quadrature_rule(&ir_ref_master);
						
						if(use_volume_differential) {
							master_fe->reinit(&el_master);
							master_fe_hack->reinit(&el_master);
						} else {
							master_fe->reinit(&el_master, side_index_master);
							master_fe_hack->reinit(&el_master, side_index_master);
						}
						
						//slave fe init
						slave_fe->attach_quadrature_rule(&ir_ref_slave);
						
						if(use_volume_differential) {
							slave_fe->reinit(&el_slave);
							master_fe_hack->reinit(&el_slave);
						} else {
							slave_fe->reinit(&el_slave, side_index_slave);
							slave_fe_hack->reinit(&el_slave, side_index_slave);
						}

						assert( check_positive_funs(*slave_fe) );
						assert( check_positive_funs(*master_fe) );

						
						//prepare result
						surface_assemble->parent_element_master  = index_master;
						
						surface_assemble->id_master 			 = el_master.id();
						
						surface_assemble->parent_element_slave   = index_slave;
						
						surface_assemble->id_slave 			     = el_slave.id();
						
						surface_assemble->coupling.zero();
						
						elemmat.zero();
						// massmat.zero();
						
						std::vector<bool> node_is_boundary_slave;
						std::vector<bool> node_is_boundary_master;
						
						nodes_are_boundary_hack(*slave_fe_hack, node_is_boundary_slave);
						nodes_are_boundary_hack(*master_fe_hack, node_is_boundary_master);

						assert( check_boundary_and_fun_consistent(*slave_fe, node_is_boundary_slave) );
						
						assert( check_node_is_boundary((*master_slave->active_local_elements_begin())->type(), node_is_boundary_master) );
						
						if(use_biorth) {
							assemble_trace_biorth_weights_from_space((*master_slave->active_local_elements_begin())->type(),
																	 node_is_boundary_slave,
																	 biorth_weights);
							
							mortar_assemble_weighted_biorth(*master_fe, *slave_fe, biorth_weights, elemmat);
						} else {
							mortar_assemble(*master_fe, *slave_fe, elemmat);
						}
						
						
						const libMesh::Point pp = side_master->point(0);
						const Real plane_offset = n_master.contract(pp);
						
						
						if(use_biorth) {
							mortar_normal_and_gap_assemble_weighted_biorth(
																		   *slave_fe,
																		   dim,
																		   n_slave,
																		   n_master,
																		   plane_offset,
																		   biorth_weights,
																		   surface_assemble->normals,
																		   surface_assemble->gap);
						} else {
							mortar_normal_and_gap_assemble(
														   dim,
														   *slave_fe,
														   n_slave,
														   n_master,
														   plane_offset,
														   surface_assemble->normals,
														   surface_assemble->gap);
						}
						//////////////////////////////////////////////////////////////////////////////////////
						
						rel_area_buffer.add(side_id_slave[side_index_slave], 0, surface_assemble->relative_area);
						
						const auto &dofs_master = master.dof_map();
						const auto &dofs_slave  = slave.dof_map();
						
						int n_nodes_face_slave  = side_slave->n_nodes();
						int n_nodes_face_master = side_master->n_nodes();
						
						
						//plot_lines(dim_master,  side_master->n_nodes(), &side_polygon_master.get_values()[0] , "master/" + std::to_string(side_id_master[0]));
						//plot_lines(dim_slave, side_slave->n_nodes(), &side_polygon_slave.get_values()[0] , "slave/" + std::to_string(side_id_slave[0]));
						
						std::vector<dof_id_type> face_node_id_slave(dofs_slave.size(),   -1);
						std::vector<dof_id_type> face_node_id_master(dofs_master.size(), -1);
						
						//generate face-node ids for slave
						int offset = 0;
						for(uint i = 0; i < node_is_boundary_slave.size(); ++i){
							if (node_is_boundary_slave[i]) {
								face_node_id_slave[i] =  side_id_slave[side_index_slave] * n_nodes_face_slave + offset++;
							}
						}
						
						//generate face-node ids for master
						offset = 0;
						for(uint i = 0; i < node_is_boundary_master.size(); ++i) {
							if (node_is_boundary_master[i]) {
								face_node_id_master[i] = side_id_master[side_index_master] * n_nodes_face_master + offset++;
							}
						}
						
						//fill-up slave permutation
						for(int i = 0; i <  dofs_slave.size(); ++i) {
							const long dof_I = dofs_slave[i];
							const long dof_J = face_node_id_slave[i];
							
							if(node_is_boundary_slave[i]) {
								P_buffer.set(dof_I, dof_J, 1.);
							}
						}
						
						//fill-up master permutation
						for(int i = 0; i <  face_node_id_master.size(); ++i) {
							const long dof_I = dofs_master[i];
							const long dof_J = face_node_id_master[i];
							
							if(node_is_boundary_master[i]) {
								Q_buffer.set(dof_I, dof_J, 1.);
							}
						}
						
						// if(enable_vis) {
						// 	DenseVector<Real> v(elemmat.n());
						// 	std::fill(v.get_values().begin(), v.get_values().end(), 1.0);
						// 	DenseVector<Real> r(elemmat.m());
						// 	std::fill(r.get_values().begin(), r.get_values().end(), 0.0);
						// 	elemmat.vector_mult(r, v);
						// 	r.print(std::cout);
						// }
						
						auto partial_sum = std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(), libMesh::Real(0.0));
						assert(partial_sum > 0);
						// assert( std::abs(surface_assemble->relative_area - 0.5*std::accumulate(ir_ref_slave.get_weights().begin(), ir_ref_slave.get_weights().end(), 0.0)) < 1e-8);

						// std::cout << "====================================\n";
						// std::cout << "areas: (" << side_id_slave[side_index_slave] << ", " << side_id_master[side_index_master] << ") "  << partial_sum << " == " << surface_assemble->isect_area << " ratio: " <<  (surface_assemble->isect_area/partial_sum) << std::endl;
						// std::cout << "====================================\n";
						assert( std::abs(partial_sum - surface_assemble->isect_area) < 1e-8);


						
						local_element_matrices_sum += partial_sum;
						
						assert(dofs_slave.size() == elemmat.m());
						assert(dofs_master.size() == elemmat.n());
						
						for(int i = 0; i < face_node_id_slave.size(); ++i) {
							if(!node_is_boundary_slave[i]) continue;
							
							const long dof_I = face_node_id_slave[i];
							
							gap_buffer.add(dof_I, 0, surface_assemble->gap(i));
							
							for (int k = 0; k < dim_slave; ++k) {
								normal_buffer.add(dof_I, k, surface_assemble->normals(i, k));
							}
							
							for(int j = 0; j <  face_node_id_master.size(); ++j) {
								if(!node_is_boundary_master[j]) continue;
								
								const long dof_J = face_node_id_master[j];
								
								B_buffer.add(dof_I, dof_J, elemmat(i, j));
							}
						}
						
						intersected = true;
					}
				}
			}
			
			return intersected;
			
		};
		
		if(!SurfaceAssemble<Dimensions>(comm, master_slave, dof_map, var_num, fun, settings, search_radius, tags, use_biorth)) {
			return false;
		}
		
		
		double volumes[1] = { local_element_matrices_sum };
		
		comm.all_reduce(volumes, 1, moonolith::MPISum());
		
		const processor_id_type master_proc_id  = master_slave->processor_id();
		
		const dof_id_type n_dofs_on_proc_master = dof_map->n_dofs_on_processor(master_proc_id);
		
		const processor_id_type slave_proc_id   = master_slave->processor_id();
		
		const dof_id_type n_dofs_on_proc_slave  = dof_map->n_dofs_on_processor(slave_proc_id);
		
		if(comm.is_root()) {
			std::cout << "sum(B_tilde): " << volumes[0] <<std::endl;
		}
		
		std::vector<moonolith::Integer>  ownership_ranges_master(comm.size()+1, 0);
		std::vector<moonolith::Integer>  ownership_ranges_slave(comm.size()+1, 0);

		
		const int n_nodes_x_face = master_slave->elem(0)->build_side_ptr(0)->n_nodes();
		std::vector<moonolith::Integer>  side_node_ownership_ranges = local_fun_spaces_new->ownershipRangesFaceID();
		
		for(SizeType i = 0; i < side_node_ownership_ranges.size(); ++i) {
			side_node_ownership_ranges[i] *= n_nodes_x_face;
		}
		
		ownership_ranges_master[comm.rank() + 1] += static_cast<unsigned int>(n_dofs_on_proc_master);
		ownership_ranges_slave[comm.rank()  + 1] += static_cast<unsigned int>(n_dofs_on_proc_slave);
		
		comm.all_reduce(&ownership_ranges_master[0], ownership_ranges_master.size(), moonolith::MPISum());
		comm.all_reduce(&ownership_ranges_slave[0],  ownership_ranges_slave.size(),  moonolith::MPISum());
		
		std::partial_sum(ownership_ranges_master.begin(), ownership_ranges_master.end(),
						 ownership_ranges_master.begin());
		
		std::partial_sum(ownership_ranges_slave.begin(), ownership_ranges_slave.end(),
						 ownership_ranges_slave.begin());
		
		
		const SizeType n_local_dofs_slave     = ownership_ranges_slave [comm.rank() + 1] - ownership_ranges_slave [comm.rank()];
		const SizeType n_local_dofs_master    = ownership_ranges_master[comm.rank() + 1] - ownership_ranges_master[comm.rank()];
		const SizeType n_local_side_node_dofs = side_node_ownership_ranges[comm.rank() + 1] - side_node_ownership_ranges[comm.rank()];
		
		B_buffer.finalize_local_structure();
		P_buffer.finalize_local_structure();
		Q_buffer.finalize_local_structure();
		rel_area_buffer.finalize_local_structure();
		gap_buffer.finalize_local_structure();
		normal_buffer.finalize_local_structure();
		
		SizeType sizes[6] = {
			B_buffer.rows(), B_buffer.cols(),
			P_buffer.rows(), P_buffer.cols(),
			Q_buffer.rows(), Q_buffer.cols()
		};
		
		comm.all_reduce(sizes, 6, moonolith::MPIMax());
		
		const SizeType n_side_node_dofs = std::max(std::max(sizes[0], sizes[1]), std::max(sizes[3], sizes[5]));
		std::cout << "n_side_node_dofs = " << n_side_node_dofs << std::endl;
		
		B_buffer.set_size(n_side_node_dofs, n_side_node_dofs);
		P_buffer.set_size(sizes[2], n_side_node_dofs);
		Q_buffer.set_size(sizes[4], n_side_node_dofs);
		rel_area_buffer.set_size(n_side_node_dofs, 1);
		gap_buffer.set_size(n_side_node_dofs, 1);
		normal_buffer.set_size(n_side_node_dofs, dim);
		
		moonolith::Redistribute< moonolith::SparseMatrix<double> > redist(comm.get_mpi_comm());
		redist.apply(side_node_ownership_ranges, B_buffer,     	moonolith::AddAssign<double>());
		redist.apply(side_node_ownership_ranges, gap_buffer,    moonolith::AddAssign<double>());
		redist.apply(side_node_ownership_ranges, normal_buffer, moonolith::AddAssign<double>());
		
		redist.apply(local_fun_spaces_new->ownershipRangesFaceID(), rel_area_buffer, moonolith::AddAssign<double>());
		
		redist.apply(ownership_ranges_slave, P_buffer,  moonolith::Assign<double>());
		redist.apply(ownership_ranges_master, Q_buffer, moonolith::Assign<double>());
		
		std::cout << "n_local_side_node_dofs: " << n_local_side_node_dofs << std::endl;
		
		moonolith::root_describe("petsc rel_area_buffer assembly begin", comm, std::cout);
		
		std::vector<double> detect_negative(n_local_side_node_dofs);
		std::vector<bool> remove_row(n_local_side_node_dofs);
		
		if(!remove_row.empty()) {
			long n_remove_rows = 0;

			std::fill(detect_negative.begin(), detect_negative.end(), 0.);
			std::fill(remove_row.begin(), remove_row.end(), false);
			
			
			//hack this is curing some symptoms of a bug but not the cause
			{
				for(auto it = B_buffer.iter(); it; ++it) {
					const SizeType index  = it.row() - side_node_ownership_ranges[comm.rank()];
					detect_negative[index] += *it;
				}
			}
			
			
			{
				for (auto it = rel_area_buffer.iter(); it; ++it) {
					bool must_remove = *it < 1 - 1e-8;
					
					if(!must_remove) {
						const SizeType faceId = it.row();
						
						for(int k = 0; k < n_nodes_x_face; ++k) {
							const SizeType nodeId = faceId * n_nodes_x_face + k;
							const SizeType index  = nodeId - side_node_ownership_ranges[comm.rank()];
							
							if(detect_negative[index] < 0.) {
								must_remove = true;
								std::cerr << "[Warning] removing element with negative contribution face id: " << faceId << ", node id: " << nodeId << ", node offset: " << k  << " value: " << detect_negative[index] << std::endl;
								break;
							}
						}
					}
					
					
					if(must_remove) {
						const SizeType faceId = it.row();
						
						for(int k = 0; k < n_nodes_x_face; ++k) {
							const SizeType nodeId = faceId * n_nodes_x_face + k;
							const SizeType index  = nodeId - side_node_ownership_ranges[comm.rank()];
							assert(index < remove_row.size());
							remove_row[index] = true;
							++n_remove_rows;
						}
					}
				}
			}
			
			// std::cout << "n_remove_rows: " <<n_remove_rows << std::endl;
		}
		
		moonolith::root_describe("petsc b_buffer assembly begin", comm, std::cout);
		
		DVectord is_contact_node_tilde = local_zeros(n_local_side_node_dofs);
		DVectord gap_tilde = local_zeros(n_local_side_node_dofs);
		{
			Write<DVectord> w_g(gap_tilde);
			Write<DVectord> w_i(is_contact_node_tilde);

			for (auto it = gap_buffer.iter(); it; ++it) {
				
				const SizeType index = it.row() - side_node_ownership_ranges[comm.rank()];
				assert(index < remove_row.size());
				
				if(!remove_row[index]) {
					gap_tilde.set(it.row(), *it);
					is_contact_node_tilde.set(it.row(), 1.0);
				}
			}
		}
		
		DSMatrixd normal_tilde = utopia::local_sparse(n_local_side_node_dofs, dim, dim);
		{
			utopia::Write<utopia::DSMatrixd> write(normal_tilde);
			for (auto it = normal_buffer.iter(); it; ++it) {
				
				const SizeType index = it.row() - side_node_ownership_ranges[comm.rank()];
				assert(index < remove_row.size());
				
				if(!remove_row[index]) {
					normal_tilde.set(it.row(), it.col(), *it);
				}
			}
		}
		
		SizeType n_max_row_entries_bpq[3] = { B_buffer.local_max_entries_x_col(), P_buffer.local_max_entries_x_col(), Q_buffer.local_max_entries_x_col() };
		comm.all_reduce(n_max_row_entries_bpq, 3, moonolith::MPIMax());
		
		const SizeType n_max_row_entries_b = n_max_row_entries_bpq[0];
		const SizeType n_max_row_entries_p = n_max_row_entries_bpq[1];
		const SizeType n_max_row_entries_q = n_max_row_entries_bpq[2];
		
		DSMatrixd B_tilde = utopia::local_sparse(n_local_side_node_dofs, n_local_side_node_dofs, n_max_row_entries_b);
		{
			utopia::Write<utopia::DSMatrixd> write(B_tilde);
			for (auto it = B_buffer.iter(); it; ++it) {
				
				const SizeType index = it.row() - side_node_ownership_ranges[comm.rank()];
				assert(index < remove_row.size());
				
				if(!remove_row[index]) {
					B_tilde.set(it.row(), it.col(), *it);
				}
			}
		}
		
		// rel_area_buffer.save("area.txt");
		// B_buffer.save("B_tilde.txt");
		
#ifndef NEBUG
		DVectord B_tilde_lumped = sum(B_tilde, 1);
		each_read(B_tilde_lumped, [](const SizeType i, const double value) {
			//			assert(value >= 0.);
		});
#endif //NEBUG
		
		// comm.barrier();
		// moonolith::root_describe("petsc p_buffer assembly begin", comm, std::cout);
		
		DSMatrixd P = utopia::local_sparse(n_local_dofs_slave, n_local_side_node_dofs, n_max_row_entries_p);
		{
			utopia::Write<utopia::DSMatrixd> write(P);
			for (auto it = P_buffer.iter(); it; ++it) {
				P.set(it.row(), it.col(), *it);
			}
		}
		
		// comm.barrier();
		// moonolith::root_describe("petsc q_buffer assembly begin", comm, std::cout);
		
		DSMatrixd Q = utopia::local_sparse(n_local_dofs_master, n_local_side_node_dofs, n_max_row_entries_q);
		{
			utopia::Write<utopia::DSMatrixd> write(Q);
			for (auto it = Q_buffer.iter(); it; ++it) {
				Q.set(it.row(), it.col(), *it);
			}
		}
		
		DSMatrixd Q_t = transpose(Q);
		DSMatrixd B_x = P * B_tilde * Q_t;
		
		normals = P * normal_tilde;
		
		DVectord gap_x = P * gap_tilde;
		
		is_contact_node = P * is_contact_node_tilde;
		
		
		DVectord normals_vec = local_zeros(n_local_dofs_slave);
		{
			Write<DVectord> w(normals_vec);
			
			each_read(normals, [&](const SizeType i, const SizeType j, const double value){
				normals_vec.set(i + j, value);
			});
		}
		
		
		bool has_contact = false;
		
		{
			Write<DVectord> w_i(is_contact_node);
			each_read(is_contact_node, [&](const SizeType i , const double value){
				if (value > 0)
				{
					is_contact_node.set(i, 1);
					has_contact = true;
				}
			});
		}
		
		orthogonal_trafos = local_sparse(n_local_dofs_slave , n_local_dofs_slave , dim);
		{
			typedef Intersector::Scalar Scalar;
			std::vector<Scalar> normal(dim, 0);
			std::vector<Scalar> H(dim * dim, 0);
			
			Read<DVectord>  r_n(normals_vec);
			Write<DSMatrixd> w_o(orthogonal_trafos);
			Read<DVectord> r_icn(is_contact_node);
			
			bool check_has_contact = false;
			
			utopia::Range r = utopia::range(normals_vec);
			for(uint i = r.begin(); i < r.end(); i += dim) {
				bool use_identity = true;
				bool is_cn_i = is_contact_node.get(i) > 0;
				
				if(is_cn_i) {
					check_has_contact = true;
					
					for(uint d = 0; d < dim; ++d) {
						normal[d] = normals_vec.get(i + d);
					}
					
					normalize(normal);
					
					if(std::abs(normal[0] - 1.) > 1e-8) {
						use_identity = false;
						
						//-e1 basis vector
						normal[0] -= 1;
						normalize(normal);
						
						if(dim == 2) {
							isector.householder_reflection_2(&normal[0], &H[0]);
						} else {
							isector.householder_reflection_3(&normal[0], &H[0]);
						}
						
						for(uint di = 0; di < dim; ++di) {
							for(uint dj = 0; dj < dim; ++dj) {
								orthogonal_trafos.set((i + di), (i + dj), H[di * dim + dj]);
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
			
			if(!check_has_contact == has_contact) {
				std::cerr << "inconsistent contact determination" << std::endl;
			}
		}
		
		auto s_gap = local_size(gap_x);
		gap = local_zeros(s_gap);
		
		static const double LARGE_VALUE = 10000;
		{
			Write<DVectord> w_g(gap);
			Read<DVectord> r_icn(is_contact_node);
			
			each_read(gap_x, [&](const SizeType i, const double value) {
				const SizeType offset = i;
				
				if(is_contact_node.get(i) > 0) {
					gap.set(offset, value);
				} else {
					gap.set(offset, LARGE_VALUE);
				}
			});
		}
		
		auto size_B_x = local_size(B_x);
		B = local_sparse(size_B_x.get(0), size_B_x.get(1), n_max_row_entries_b * dim);
		
		{
			Write<DSMatrixd> w_B(B);
			each_read(B_x, [&](const SizeType i, const SizeType j, const double value) {
				for(SizeType d = 0; d < dim; ++d) {
					B.set(i + d, j + d, value);
				}
			});
		}
		
		// disp("B_tilde:");
		// disp(B.size());
		
		// disp("Q:");
		// disp(Q_t.size());
		
		// disp("P:");
		// disp(P.size());
		
		// disp("B:");
		// disp(B.size());
		
		
		// write("B_tilde.m", B_tilde);
		// write("Q.m", Q);
		// write("P.m", P);
		
		// write("B.m", B);
		// write("O.m", orthogonal_trafos);
		// write("g.m", gap);
		// write("c.m", is_contact_node);
		
		comm.barrier();
		moonolith::root_describe("contact assembly end", comm, std::cout);
		return true;
	}
	
	bool assemble_contact(moonolith::Communicator &comm,
						  const std::shared_ptr<MeshBase> &mesh,
						  const std::shared_ptr<DofMap> &dof_map,
						  const unsigned int var_num,
						  DSMatrixd &B,
						  DSMatrixd &orthogonal_trafos,
						  DVectord &gap,
						  DSMatrixd &normals,
						  DVectord &is_contact_node,
						  const libMesh::Real search_radius,
						  const std::vector< std::pair<int, int> > &tags,
						  const bool use_biorth,
						  const bool use_volume_differential)
	{
		moonolith::SearchSettings settings;
		// settings.verbosity_level = 3;

		if(mesh->mesh_dimension() == 2) {
			return utopia::SurfaceAssemble<2>(comm, mesh, dof_map, var_num, B,  orthogonal_trafos, gap, normals, is_contact_node, settings, search_radius, tags, use_biorth, use_volume_differential);
		}
		
		
		if(mesh->mesh_dimension() == 3) {
			return utopia::SurfaceAssemble<3>(comm, mesh, dof_map, var_num, B,  orthogonal_trafos, gap, normals,  is_contact_node, settings, search_radius, tags, use_biorth, use_volume_differential);
		}
		
		assert(false && "Dimension not supported!");
		return false;
	}
	
	bool assemble_contact(moonolith::Communicator &comm,
						  const std::shared_ptr<MeshBase> &mesh,
						  const std::shared_ptr<DofMap> &dof_map,
						  const unsigned int var_num,
						  DSMatrixd &B,
						  DSMatrixd &orthogonal_trafos,
						  DVectord &gap,
						  DSMatrixd &normals,
						  DVectord &is_contact_node,
						  const libMesh::Real search_radius,
						  const int tag_1,
						  const int tag_2,
						  const bool use_biorth,
						  const bool use_volume_differential)
	{
		return assemble_contact(
								comm, mesh, dof_map, var_num,
								B, orthogonal_trafos,
								gap, normals, is_contact_node,
								search_radius,
								{ {tag_1, tag_2} },
								use_biorth,
								use_volume_differential);
	}
	
	void convert_normal_matrix_to_vector(const DSMatrixd &mat, DVectord &vec)
	{
		auto s = local_size(mat);
		
		vec = local_zeros(s.get(0));
		DVectord norms = local_zeros(s.get(0));
		
		auto s_ns = local_size(norms);
		
		{
			Write<DVectord> w_ns(norms);
			each_read(mat, [&](const SizeType i, const SizeType j, const double value){
				norms.add(i, value * value);
			});
		}
		
		norms = sqrt(norms);
		
		{
			Write<DVectord> w(vec);
			Read<DVectord> r_ns(norms);
			
			each_read(mat, [&](const SizeType i, const SizeType j, const double value){
				vec.set(i + j, value/norms.get(i));
			});
		}
	}
	
	bool assemble_contact(
						  moonolith::Communicator &comm,
						  const std::shared_ptr<libMesh::MeshBase> &mesh,
						  const std::shared_ptr<libMesh::DofMap> &dof_map,
						  const unsigned int var_num,
						  DSMatrixd &B,
						  DSMatrixd &orthogonal_trafos,
						  DVectord &gap,
						  DVectord &normals,
						  DVectord &is_contact_node,
						  const libMesh::Real search_radius,
						  const std::vector< std::pair<int, int> > &tags,
						  const bool use_biorth,
						  const bool use_volume_differential)
	{
		
		DSMatrixd direction_matrix;
		if(!assemble_contact(
							 comm, mesh, dof_map, var_num,
							 B, orthogonal_trafos,
							 gap, direction_matrix, is_contact_node,
							 search_radius,
							 tags,
							 use_biorth,
							 use_volume_differential)) {
			return false;
		}
		
		convert_normal_matrix_to_vector(direction_matrix, normals);
		return true;
	}
}

#include "utopia_assemble_contact.hpp"

#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/serial_mesh.h"

#include "express_Profiler.hpp"
#include "express_Redistribute.hpp"

#include "utopia_fe.hpp"
#include "MortarAssemble.hpp"

#include "MortarAssemble.hpp"
#include "MortarAssembler.hpp"

#include "utopia_Socket.hpp"
#include "utopia_STree.hpp"
#include "utopia_ElementDofMap.hpp"
#include "utopia_FESpaceAdapter.hpp"

#include "MapSparseMatrix.hpp"

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
	
	//Hacky but it is the only way with libmesh that we know of without having to clone the boundary info!!!!
	inline static void nodes_are_boundary_hack(
											   const libMesh::DenseMatrix<libMesh::Real> &mat,
											   std::vector<bool> &rows,
											   std::vector<bool> &cols)
	{
		
		rows.resize(mat.m());
		cols.resize(mat.n());
		
		std::fill(rows.begin(), rows.end(), 0);
		std::fill(cols.begin(), cols.end(), 0);
		
		std::vector<double> sum_rows(mat.m(), 0.);
		std::vector<double> sum_cols(mat.n(), 0.);
		double sum_all = 0.0;
		double vol_check = 0.0;
		
#ifndef NDEBUG
		std::vector<double> sum_rows_with_sign(mat.m(), 0.);
		
#endif //NDEBUG
		
		for(int i = 0; i < sum_rows.size(); ++i) {
			for(int j = 0; j < sum_cols.size(); ++j) {
				const double abs_ij = std::abs(mat(i, j));
				sum_rows[i] += abs_ij;
				sum_cols[j] += abs_ij;
				sum_all += abs_ij;
				
				vol_check += mat(i, j);
				
#ifndef NDEBUG
				sum_rows_with_sign[i] += mat(i, j);
#endif //NDEBUG
			}
		}
		
		for(int i = 0; i < sum_rows.size(); ++i) {
			if(std::abs(sum_rows[i]/sum_all) > 1e-8) {
				rows[i] = true;
			}
		}
		
		for(int i = 0; i < sum_cols.size(); ++i) {
			if(std::abs(sum_cols[i]/sum_all) > 1e-8) {
				cols[i] = true;
			}
		}
		
		assert(vol_check > 0);
#ifndef NDEBUG
		for(int i = 0; i < sum_rows_with_sign.size(); ++i) {
			assert(sum_rows_with_sign[i]/sum_all >= -1e-8);
		}
#endif //NDEBUG
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
							cutk::OutputStream &os)
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
		os.requestSpace( (n_elements * 8 + n_nodes * dim) * (sizeof(double) + sizeof(long)) );
		
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
			
			const int k = *it;
			
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
			
			bool  size=true;
			
			int volume_tag;
			
			volume_tag=subdomain_id[elem->id()].global.at(0);
			
			os << volume_tag;
			
			int side_set_tag;
			
			int face_id;
			
			bool check_side_id_one = true;
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
										cutk::OutputStream &os)
	{
		write_space(begin, end, *fespace.mesh(), fespace.dof_map(), fespace.variable_number(), fespace.variable_number(), fespace.subdomain_id(), fespace.side_set_id(), fespace.face_set_id_global(), os);
	}
	
	
	static void read_space(cutk::InputStream &is, cutk::shared_ptr<MeshBase> & space,
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
			
			int volume_tag, side_set_tag, face_id;
			
			bool on_boundary=false;
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
	
	static void read_spaces(cutk::InputStream &is, FESpaceAdapter &utopiamesh, const libMesh::Parallel::Communicator &comm_mesh)
	{
		read_space(is, utopiamesh.mesh(), utopiamesh.dof_map(), utopiamesh.variable_number(), utopiamesh.variable_order(), utopiamesh.subdomain_id(), utopiamesh.side_set_id(), utopiamesh.face_set_id_global(), comm_mesh);
	}
	
	
	template<int Dimensions, class Fun>
	static bool SurfaceAssemble(express::Communicator &comm,
								const std::shared_ptr<MeshBase> &master_slave,
								const std::shared_ptr<DofMap> &dof_map,
								const unsigned int var_num,
								Fun process_fun,
								const cutk::Settings &settings,
								const libMesh::Real search_radius,
								const std::vector< std::pair<int, int> >  &tags,
								const bool use_biorth)
	{
		
		std::shared_ptr<FESpaceAdapter> local_fun_spaces_new = cutk::make_shared<FESpaceAdapter>(master_slave, dof_map, var_num, tags);
		auto predicate = std::make_shared<cutlibpp::MasterAndSlave>();
		
		for(auto t : tags)
			predicate->add(t.first, t.second);
		
		using namespace cutlibpp;
		using namespace express;
		using namespace cutk;
		
		typedef STree<Dimensions> NTreeT;
		typedef typename NTreeT::DataContainer DataContainer;
		typedef typename NTreeT::DataType SurfaceAdapter;
		
		long maxNElements = 40;
		long maxDepth = 5;
		
		
		if (!settings.get("max_depth").isNull()) {
			maxDepth = settings.get("max_depth").toInt();
		}
		
		const auto &mesh = master_slave;
		
		const int n_elements = mesh->n_elem();
		
		const Parallel::Communicator &libmesh_comm_mesh = master_slave->comm();
		
		const int dim_master = master_slave->mesh_dimension();
		const int dim_slave = master_slave->mesh_dimension();
		
		MeshBase::const_element_iterator e_it = mesh->active_elements_begin();
		const MeshBase::const_element_iterator e_end = mesh->active_elements_end();
		std::vector<int> block_id;
		std::vector<int> block_id_def;
		
		EXPRESS_EVENT_BEGIN("create_adapters");
		////////////////////////////////////////////////////////////////////////////////////////////////////
		cutk::shared_ptr<NTreeT> tree = NTreeT::New(predicate, maxNElements, maxDepth);
		tree->reserve(n_elements);
		
		std::cout << "nElements = tree->memory().nData()_inside " << n_elements << std::endl;
		
		std::shared_ptr<FESpaceAdapter> local_spaces = make_shared<FESpaceAdapter>(master_slave, dof_map, var_num, tags);
		
		int jj=0;
		
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
		
		tree->getRoot()->getBound().staticBound().enlarge(1e-8);
		
		
		////////////////////////////////////////////////////////////////////////////////////////////////////
		EXPRESS_EVENT_END("create_adapters");
		
		//Just to have an indexed-storage
		std::map<long, cutk::shared_ptr<FESpaceAdapter> > utopiamesh;
		std::map<long, std::vector<cutk::shared_ptr<FESpaceAdapter> > > migrated_meshes;
		
		
		auto read = [&utopiamesh, &migrated_meshes, block_id, comm, &libmesh_comm_mesh, search_radius]
		(
		 const long ownerrank,
		 const long senderrank,
		 bool is_forwarding, DataContainer &data,
		 InputStream &in
		 ) {
			
			CHECK_STREAM_READ_BEGIN("vol_proj", in);
			
			cutk::shared_ptr<FESpaceAdapter> proc_space = cutk::make_shared<FESpaceAdapter>(comm);
			
			read_spaces(in, *proc_space, libmesh_comm_mesh);
			
			if (!is_forwarding) {
				assert(!utopiamesh[ownerrank]);
				utopiamesh[ownerrank] = proc_space;
			} else {
				migrated_meshes[ownerrank].push_back(proc_space);
			}
			
			data.reserve(data.size() + proc_space->n_elements());
			
			auto s = proc_space->mesh();
			
			int i=0;
			for (int i = 0; i<s->n_elem(); ++i) {
				auto elem=s->elem(i);
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
				cutk::shared_ptr<FESpaceAdapter> spaceptr = it->second;
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
		
		
		cutk::Settings custom_settings = settings;
		// custom_settings.set("disable_redistribution", cutk::Boolean(true));
		custom_settings.set("verbosity_level", cutk::Integer(1));
		
		cutlibpp::search_and_compute(comm, tree, predicate, read, write, fun, custom_settings);
		
		long n_total_candidates = n_projections + n_false_positives;
		
		long n_collection[3] = {n_projections, n_total_candidates, n_false_positives};
		comm.allReduce(n_collection, 3, express::MPISum());
		
		if (comm.isRoot()) {
			std::cout << "n_intersections: " << n_collection[0]
			<< ", n_total_candidates: " 	 << n_collection[1]
			<< ", n_false_positives: " 	     << n_collection[2] << std::endl;
		}
		
		return true;
	}
	
	template<int Dimensions>
	bool SurfaceAssemble(
						 express::Communicator &comm,
						 const std::shared_ptr<MeshBase> &master_slave,
						 const std::shared_ptr<DofMap> &dof_map,
						 const unsigned int var_num,
						 DSMatrixd &B,
						 DSMatrixd &orthogonal_trafos,
						 DVectord &gap,
						 DSMatrixd &normals,
						 DVectord &is_contact_node,
						 const cutk::Settings &settings,
						 const libMesh::Real search_radius,
						 const std::vector< std::pair<int, int> > &tags,
						 const bool use_biorth)
	{
		std::shared_ptr<FESpaceAdapter> local_fun_spaces_new = cutk::make_shared<FESpaceAdapter>(master_slave, dof_map, var_num, tags);
		
		libMesh::DenseMatrix<libMesh::Real> src_pts;
		libMesh::DenseMatrix<libMesh::Real> dest_pts;
		libMesh::DenseMatrix<libMesh::Real> intersection2;
		Polyhedron src_poly, dest_poly;
		Polyhedron  intersection3,temp_poly;
		Intersector isector;
		
		std::shared_ptr<MeshBase> master_slave_space = master_slave;
		
		auto predicate = std::make_shared<cutlibpp::MasterAndSlave>();
		
		for(auto t : tags){
			predicate->add(t.first, t.second);
			std::cout << "[Status] Added master-slave pair = "<< t.first << ", " << t.second << std::endl;
		}
		
		static const double tol = 1e-8;
		
		std::vector<libMesh::dof_id_type> master_dofs, slave_dofs;
		
		libMesh::DenseMatrix<libMesh::Real> elemmat;
		libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;
		
		DenseMatrix<Real> side_polygon_master, side_polygon_slave;
		DenseMatrix<Real> isect_polygon_master, isect_polygon_slave;
		
		std::shared_ptr<Transform> trafo_master;
		std::shared_ptr<Transform> trafo_slave;
		
		Point n_master, n_slave;
		
		const int dim = master_slave->mesh_dimension();
		
		libMesh::Real total_intersection_volume = 0.0;
		libMesh::Real local_element_matrices_sum = 0.0;
		
		//all face dof buffers
		express::MapSparseMatrix<double> B_buffer;
		express::MapSparseMatrix<double> P_buffer;
		express::MapSparseMatrix<double> Q_buffer;
		express::MapSparseMatrix<double> rel_area_buffer;
		express::MapSparseMatrix<double> gap_buffer;
		express::MapSparseMatrix<double> normal_buffer;
		
		// std::cout<<"*********** master_slave->dof_map().n_dofs() = "<<  dof_map->n_dofs() <<std::endl;
		
		DenseMatrix<Real> biorth_weights;
		
		auto fun = [&](const SElementAdapter<Dimensions> &master,
					   const SElementAdapter<Dimensions> &slave) -> bool {
			
			long n_intersections = 0;
			
			using namespace cutlibpp;
			using namespace express;
			using namespace cutk;
			
			const auto &master_mesh = master.space();
			const auto &slave_mesh  = slave.space();
			
			libMesh::DenseMatrix<libMesh::Real> elemmat;
			
			const int src_index  = master.element();
			const int dest_index = slave.element();
			
			auto &el_master  = *master_mesh.elem(src_index);
			auto &dest_el = *slave_mesh.elem(dest_index);
			
			const int dim_master = master_mesh.mesh_dimension();
			const int dim_slave = slave_mesh.mesh_dimension();
			
			Box box_master(dim_master), box_slave(dim_slave);
			
			QMortar src_ir_ref(dim_master);
			QMortar src_ir(dim_master);
			QMortar dest_ir(dim_slave);
			QMortar dest_ir_ref(dim_slave);
			
			//only works because there are not mixed elements
			const int approx_order = local_fun_spaces_new->variable_order()[0];
			
			std::shared_ptr<Contact> surface_assemble;
			
			const auto &side_id_master = master.dof_map_face();
			const auto &face_id_slave  = slave.dof_map_face();
			
			std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
			
			master_fe = libMesh::FEBase::build(master_mesh.mesh_dimension(), FIRST);
			slave_fe  = libMesh::FEBase::build(slave_mesh.mesh_dimension(),  FIRST);
			
			typedef Intersector::Scalar Scalar;
			
			if(dim_slave == 2)  {
				make_polygon(el_master,   src_pts);
				make_polygon(dest_el,  dest_pts);
				trafo_master  = std::make_shared<AffineTransform2>(el_master);
				trafo_slave = std::make_shared<AffineTransform2>(dest_el);
				
			}
			
			else if(dim_slave == 3) {
				make_polyhedron(el_master,  src_poly);
				make_polyhedron(dest_el, dest_poly);
				trafo_master  = std::make_shared<AffineTransform3>(el_master);
				trafo_slave = std::make_shared<AffineTransform3>(dest_el);
				
			}
			
			bool intersected = false;
			
			for(uint side_index_master = 0;
				side_index_master < el_master.n_sides();
				++side_index_master) {
				
				if(side_id_master[side_index_master] < 0) continue;
				
				if(el_master.neighbor_ptr(side_index_master) != nullptr) continue;
				
				auto side_master = el_master.build_side_ptr(side_index_master);
				
				compute_side_normal(dim_master, *side_master, n_master);
				fix_normal_orientation(el_master, side_index_master, n_master);
				
				box_master.reset();
				enlarge_box_from_side(dim_master, *side_master, box_master, search_radius);
				
				if(dim_slave == 2) {
					make_polygon(*side_master, side_polygon_master);
				} else if(dim_master == 3) {
					make_polygon_3(*side_master, side_polygon_master);
				} else {
					assert(false);
				}
				
				for(uint side_2 = 0; side_2 < dest_el.n_sides(); ++side_2) {
					if(face_id_slave[side_2] < 0) continue;
					if(dest_el.neighbor_ptr(side_2) != nullptr) continue;
					// if (!predicate->tagsAreRelated(tag_1, tag_2)) continue;
					
					auto side_slave = dest_el.build_side_ptr(side_2);
					
					compute_side_normal(dim_slave, *side_slave, n_slave);
					
					fix_normal_orientation(dest_el, side_2, n_slave);
					
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
						//plot_lines(2, 2, &side_polygon_slave.get_values()[0], "in_slave/" + std::to_string(face_id_slave[0]) + "_" + std::to_string(cos_angle));
						
						if(!project_2D(side_polygon_master, side_polygon_slave, isect_polygon_master, isect_polygon_slave)){
							continue;
						}
						
						
						const Scalar dx = side_polygon_slave(0, 0) - side_polygon_slave(1, 0);
						const Scalar dy = side_polygon_slave(0, 1) - side_polygon_slave(1, 1);
						
						const Scalar isect_dx = isect_polygon_slave(0, 0) - isect_polygon_slave(1, 0);
						const Scalar isect_dy = isect_polygon_slave(0, 1) - isect_polygon_slave(1, 1);
						
						const Scalar area   = std::sqrt(isect_dx*isect_dx + isect_dy*isect_dy);
						const Scalar weight = area/std::sqrt(dx*dx + dy*dy);
						
						if(weight < 1e-15) continue;
						
						const int order = order_for_l2_integral(dim_master, el_master, approx_order, dest_el, approx_order);
						
						make_composite_quadrature_on_surf_2D(isect_polygon_master, weight, order, src_ir);
						
						make_composite_quadrature_on_surf_2D(isect_polygon_slave, weight, order, dest_ir);
						
						pair_intersected = true;
						
						surface_assemble = std::make_shared<Contact>();
						surface_assemble->isect_area	   = area;
						surface_assemble->relative_area    = weight;
						
						// plot_polygon(2, 2, &side_polygon_master.get_values()[0], "master");
						// plot_polygon(2, 2, &side_polygon_slave.get_values()[0], "slave");
						
						
					} else if(dim_slave == 3) {
						make_polygon_3(*side_slave, side_polygon_slave);
						
						if(!project_3D(
									   side_polygon_master,
									   side_polygon_slave,
									   isect_polygon_master,
									   isect_polygon_slave))
						{
							continue;
						}
						
						const Scalar area_slave = isector.polygon_area_3(side_polygon_slave.m(),  &side_polygon_slave.get_values()[0]);
						const Scalar area   	= isector.polygon_area_3(isect_polygon_slave.m(), &isect_polygon_slave.get_values()[0]);
						const Scalar weight 	= area/area_slave;
						
						assert(area_slave > 0);
						assert(area > 0);
						assert(weight > 0);
						
						const int order = order_for_l2_integral(dim_master, el_master, approx_order, dest_el, approx_order);
						
						make_composite_quadrature_on_surf_3D(isect_polygon_master, weight, order, src_ir);
						make_composite_quadrature_on_surf_3D(isect_polygon_slave, weight, order, dest_ir);
						
						pair_intersected = true;
						
						surface_assemble = std::make_shared<Contact>();
						surface_assemble->isect_area	= area;
						surface_assemble->relative_area = weight;
						
					} else {
						assert(false);
						return false;
					}
					
					
					if(pair_intersected) {
						
						// plot_polygon(3, isect_polygon_master.get_values().size()/3, &isect_polygon_master.get_values()[0], "master");
						// plot_polygon(3, isect_polygon_slave.get_values().size()/3, &isect_polygon_slave.get_values()[0], "slave");
						
						// std::cout << "isect: " << master.handle() << " -> " << slave.handle() << std::endl;
						
						//////////////////////////////////ASSEMBLY ////////////////////////////////////////
						//////////////////////////////////////////////////////////////////////////////////////
						transform_to_reference_surf(*trafo_master,  el_master.type(),  src_ir, src_ir_ref);
						transform_to_reference_surf(*trafo_slave, dest_el.type(), dest_ir, dest_ir_ref);
						
						master_fe->attach_quadrature_rule(&src_ir_ref);
						
						master_fe->get_phi();
						master_fe->reinit(&el_master);
						
						slave_fe->attach_quadrature_rule(&dest_ir_ref);
						
						slave_fe->get_xyz();
						slave_fe->reinit(&dest_el);
						
						
						surface_assemble->parent_element_master  = src_index;
						
						surface_assemble->id_master 			 = el_master.id();
						
						surface_assemble->parent_element_slave   = dest_index;
						
						surface_assemble->id_slave 			     = dest_el.id();
						
						surface_assemble->coupling.zero();
						
						elemmat.zero();
						
						mortar_assemble(*master_fe, *slave_fe, elemmat);
						
						std::vector<bool> node_is_boundary_slave;
						std::vector<bool> node_is_boundary_master;
						
						nodes_are_boundary_hack(elemmat, node_is_boundary_slave, node_is_boundary_master);
						
						assert( check_node_is_boundary((*master_slave->active_local_elements_begin())->type(), node_is_boundary_master) );
						
						if(use_biorth) {
							assemble_trace_biorth_weights_from_space((*master_slave->active_local_elements_begin())->type(),
																	 node_is_boundary_slave,
																	 biorth_weights);
							elemmat.zero();
							mortar_assemble_weighted_biorth(*master_fe, *slave_fe, biorth_weights, elemmat);
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
						
						rel_area_buffer.add(face_id_slave[side_2], 0, surface_assemble->relative_area);
						
						const auto &master_dofs = master.dof_map();
						const auto &slave_dofs  = slave.dof_map();
						
						int n_nodes_face_slave  = side_slave->n_nodes();
						int n_nodes_face_master = side_master->n_nodes();
						
						
						//plot_lines(dim_master,  side_master->n_nodes(), &side_polygon_master.get_values()[0] , "master/" + std::to_string(side_id_master[0]));
						//plot_lines(dim_slave, side_slave->n_nodes(), &side_polygon_slave.get_values()[0] , "slave/" + std::to_string(face_id_slave[0]));
						
						std::vector<dof_id_type> face_node_id_slave(slave_dofs.size(),   -1);
						std::vector<dof_id_type> face_node_id_master(master_dofs.size(), -1);
						
						//generate face-node ids for slave
						int offset = 0;
						for(uint i = 0; i < node_is_boundary_slave.size(); ++i){
							if (node_is_boundary_slave[i]) {
								face_node_id_slave[i] =  face_id_slave[side_2] * n_nodes_face_slave + offset++;
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
						for(int i = 0; i <  slave_dofs.size(); ++i) {
							const long dof_I = slave_dofs[i];
							const long dof_J = face_node_id_slave[i];
							
							if(node_is_boundary_slave[i]) {
								P_buffer.setAt(dof_I, dof_J, 1.);
							}
						}
						
						//fill-up master permutation
						for(int i = 0; i <  face_node_id_master.size(); ++i) {
							const long dof_I = master_dofs[i];
							const long dof_J = face_node_id_master[i];
							
							if(node_is_boundary_master[i]) {
								Q_buffer.setAt(dof_I, dof_J, 1.);
							}
						}
						
						auto partial_sum = std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(), libMesh::Real(0.0));
						assert(partial_sum > 0);
						
						local_element_matrices_sum += partial_sum;
						
						assert(slave_dofs.size() == elemmat.m());
						assert(master_dofs.size() == elemmat.n());
						
						for(int i = 0; i < face_node_id_slave.size(); ++i) {
							if(!node_is_boundary_slave[i]) continue;
							
							const long dof_I = face_node_id_slave[i];
							
							gap_buffer.add(dof_I, 0, surface_assemble->gap(i));
							
							for (int k = 0; k < dim_slave; ++k) {
								normal_buffer.add(dof_I, k, surface_assemble->normals(i,k));
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
		
		comm.allReduce(volumes, 1, express::MPISum());
		
		const processor_id_type master_proc_id  = master_slave->processor_id();
		
		const dof_id_type n_dofs_on_proc_master = dof_map->n_dofs_on_processor(master_proc_id);
		
		const processor_id_type slave_proc_id   = master_slave->processor_id();
		
		const dof_id_type n_dofs_on_proc_slave  = dof_map->n_dofs_on_processor(slave_proc_id);
		
		if(comm.isRoot()) {
			std::cout << "sum(B_tilde): " << volumes[0] <<std::endl;
		}
		
		express::Array<express::SizeType>  ownership_ranges_master(comm.size()+1);
		ownership_ranges_master.allSet(0);
		
		express::Array<express::SizeType>  ownership_ranges_slave(comm.size()+1);
		ownership_ranges_slave.allSet(0);
		
		const int n_nodes_x_face = master_slave->elem(0)->build_side_ptr(0)->n_nodes();
		express::Array<express::SizeType>  side_node_ownership_ranges = local_fun_spaces_new->ownershipRangesFaceID();
		
		for(SizeType i = 0; i < side_node_ownership_ranges.size(); ++i) {
			side_node_ownership_ranges[i] *= n_nodes_x_face;
		}
		
		ownership_ranges_master[comm.rank() + 1] += static_cast<unsigned int>(n_dofs_on_proc_master);
		ownership_ranges_slave[comm.rank()  + 1] += static_cast<unsigned int>(n_dofs_on_proc_slave);
		
		comm.allReduce(&ownership_ranges_master[0], ownership_ranges_master.size(), express::MPISum());
		comm.allReduce(&ownership_ranges_slave[0],  ownership_ranges_slave.size(),  express::MPISum());
		
		std::partial_sum(ownership_ranges_master.begin(), ownership_ranges_master.end(),
						 ownership_ranges_master.begin());
		
		std::partial_sum(ownership_ranges_slave.begin(), ownership_ranges_slave.end(),
						 ownership_ranges_slave.begin());
		
		
		const SizeType n_local_dofs_slave     = ownership_ranges_slave [comm.rank() + 1] - ownership_ranges_slave [comm.rank()];
		const SizeType n_local_dofs_master    = ownership_ranges_master[comm.rank() + 1] - ownership_ranges_master[comm.rank()];
		const SizeType n_local_side_node_dofs = side_node_ownership_ranges[comm.rank() + 1] - side_node_ownership_ranges[comm.rank()];
		
		B_buffer.finalizeStructure();
		P_buffer.finalizeStructure();
		Q_buffer.finalizeStructure();
		rel_area_buffer.finalizeStructure();
		gap_buffer.finalizeStructure();
		normal_buffer.finalizeStructure();
		
		SizeType sizes[6] = {
			B_buffer.rows(), B_buffer.columns(),
			P_buffer.rows(), P_buffer.columns(),
			Q_buffer.rows(), Q_buffer.columns()
		};
		
		comm.allReduce(sizes, 6, express::MPIMax());
		
		const SizeType n_side_node_dofs = express::Math<SizeType>::Max(sizes[0], sizes[1], sizes[3], sizes[5]);
		std::cout << "n_side_node_dofs = " << n_side_node_dofs << std::endl;
		
		B_buffer.setSize(n_side_node_dofs, n_side_node_dofs);
		P_buffer.setSize(sizes[2], n_side_node_dofs);
		Q_buffer.setSize(sizes[4], n_side_node_dofs);
		rel_area_buffer.setSize(n_side_node_dofs, 1);
		gap_buffer.setSize(n_side_node_dofs, 1);
		normal_buffer.setSize(n_side_node_dofs, dim);
		
		express::Redistribute< express::MapSparseMatrix<double> > redist(comm.getMPIComm());
		redist.apply(side_node_ownership_ranges, B_buffer,     express::AddAssign<double>());
		redist.apply(side_node_ownership_ranges, gap_buffer,       express::AddAssign<double>());
		redist.apply(side_node_ownership_ranges, normal_buffer,    express::AddAssign<double>());
		
		redist.apply(local_fun_spaces_new->ownershipRangesFaceID(), rel_area_buffer, express::AddAssign<double>());
		
		redist.apply(ownership_ranges_slave, P_buffer, express::Assign<double>());
		redist.apply(ownership_ranges_master, Q_buffer, express::Assign<double>());
		
		std::cout << "n_local_side_node_dofs: " << n_local_side_node_dofs << std::endl;
		
		express::RootDescribe("petsc rel_area_buffer assembly begin", comm, std::cout);
		
		express::Array<bool> remove_row(n_local_side_node_dofs);
		
		if(!remove_row.isNull()) {
			long n_remove_rows = 0;
			remove_row.allSet(false);
			
			{
				for (auto it = rel_area_buffer.iter(); it; ++it) {
					if(*it < 1 - 1e-8) {
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
		
		express::RootDescribe("petsc B_buffer assembly begin", comm, std::cout);
		
		DVectord is_contact_node_tilde = local_zeros(n_local_side_node_dofs);
		DVectord gap_tilde = local_zeros(n_local_side_node_dofs);
		{
			utopia::Write<utopia::DVectord> write(gap_tilde);
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
		
		SizeType n_max_row_entries_bpq[3] = { B_buffer.maxEntriesXCol(), P_buffer.maxEntriesXCol(), Q_buffer.maxEntriesXCol() };
		comm.allReduce(n_max_row_entries_bpq, 3, express::MPIMax());
		
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
		
		// comm.barrier();
		// express::RootDescribe("petsc P_buffer assembly begin", comm, std::cout);
		
		DSMatrixd P = utopia::local_sparse(n_local_dofs_slave, n_local_side_node_dofs, n_max_row_entries_p);
		{
			utopia::Write<utopia::DSMatrixd> write(P);
			for (auto it = P_buffer.iter(); it; ++it) {
				P.set(it.row(), it.col(), *it);
			}
		}
		
		// comm.barrier();
		// express::RootDescribe("petsc Q_buffer assembly begin", comm, std::cout);
		
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
		
		each_read(is_contact_node, [&](const SizeType i , const double value){
			if (value > 0)
			{
				is_contact_node.set(i, 1);
				has_contact = true;
			}
		});
		
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
		express::RootDescribe("Contact assembly end", comm, std::cout);
		return true;
	}
	
	bool assemble_contact(express::Communicator &comm,
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
						  const bool use_biorth)
	{
		cutk::Settings settings;
		if(mesh->mesh_dimension() == 2) {
			return utopia::SurfaceAssemble<2>(comm, mesh, dof_map, var_num, B,  orthogonal_trafos, gap, normals, is_contact_node, settings, search_radius, tags, use_biorth);
		}
		
		
		if(mesh->mesh_dimension() == 3) {
			return utopia::SurfaceAssemble<3>(comm, mesh, dof_map, var_num, B,  orthogonal_trafos, gap, normals,  is_contact_node, settings, search_radius, tags, use_biorth);
		}
		
		assert(false && "Dimension not supported!");
		return false;
	}
	
	bool assemble_contact(express::Communicator &comm,
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
						  const bool use_biorth)
	{
		return assemble_contact(
								comm, mesh, dof_map, var_num, 
								B, orthogonal_trafos, 
								gap, normals, is_contact_node, 
								search_radius,
								{ {tag_1, tag_2} },
								use_biorth);
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
		express::Communicator &comm,
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
		const bool use_biorth) 
	{

		DSMatrixd direction_matrix;
		if(!assemble_contact(
				comm, mesh, dof_map, var_num, 
				B, orthogonal_trafos, 
				gap, direction_matrix, is_contact_node, 
				search_radius,
				tags,
				use_biorth)) {
			return false;
		} 
		
		// auto s = local_size(gap);
		// DVectord directions = local_zeros(s.get(0));
		
		// {
		// 	Write<DVectord> w(directions);
			
		// 	each_read(direction_matrix, [&](const SizeType i, const SizeType j, const double value){
		// 		directions.set(i + j, value);
		// 	});
		// }

		// normals = local_zeros(s.get(0));

		// auto r = range(directions);
		
		// const SizeType dims = mesh->mesh_dimension();

		// std::vector<Real> n(dims);
		// {	
		// 	Read<DVectord> r_d(directions);
		// 	Read<DVectord> r_icn(is_contact_node);

		// 	for(SizeType i = r.begin(); i < r.end(); i += dims) {
		// 		if(!is_contact_node.get(i)) continue;

		// 		Real norm = 0;
		// 		for(SizeType j = 0; j < dims; ++j) {
		// 			n[j] = directions.get(i + j);
		// 			norm += n[j] * n[j];
		// 		}

		// 		norm = std::sqrt(norm);

		// 		if(norm < 1e-16) {
		// 			std::cerr << "[Warning] 0-director" << std::endl;
		// 			continue;
		// 		}

		// 		for(SizeType j = 0; j < dims; ++j) {
		// 			n[j] /= norm;
		// 			normals.set(i + j, n[j]);
		// 		}
		// 	}
		// }

		convert_normal_matrix_to_vector(direction_matrix, normals);
		
		return true;
	}
}

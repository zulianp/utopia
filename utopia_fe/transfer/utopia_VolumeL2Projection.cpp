#include "utopia_VolumeL2Projection.hpp"

#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/serial_mesh.h"


#include "Box.hpp"
#include "utopia_fe_core.hpp"
#include "MortarAssemble.hpp"
#include "utopia_BoxAdapter.hpp"
#include "utopia_VElementAdapter.hpp"
#include "utopia_VTree.hpp"
#include "utopia_ElementDofMap.hpp"
#include "utopia_FESpacesAdapter.hpp"

#include "moonolith_profiler.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_tree.hpp"
#include "moonolith_n_tree_mutator_factory.hpp"
#include "moonolith_n_tree_with_span_mutator_factory.hpp"
#include "moonolith_n_tree_with_tags_mutator_factory.hpp"
#include "moonolith_sparse_matrix.hpp"


#include <cmath>
#include <queue>
#include <algorithm>
#include <sstream>
#include <numeric>

namespace utopia {

	class FETransfer::Buffers {
	public:
		typedef libMesh::DenseMatrix<libMesh::Real> LocalMatrixT;
		typedef moonolith::SparseMatrix<double> GlobalMatrixT;

		//fem
		std::unique_ptr<libMesh::FEBase> master_fe; 
		std::unique_ptr<libMesh::FEBase> slave_fe;
		std::shared_ptr<Transform> master_trafo;
		std::shared_ptr<Transform> slave_trafo;

		//matrices
		LocalMatrixT master_pts;
		LocalMatrixT slave_pts;
		LocalMatrixT elemmat;
		LocalMatrixT cumulative_elemmat;

		std::shared_ptr<GlobalMatrixT> coupling_buffer;
		std::shared_ptr<GlobalMatrixT> mass_buffer;

		//intersections
		Intersector isector;

		Polyhedron master_poly;
		Polyhedron slave_poly;
		Polyhedron intersection3;
		LocalMatrixT intersection2;

		void init(MPI_Comm comm)
		{
			coupling_buffer = std::make_shared<GlobalMatrixT>(comm);
			mass_buffer     = std::make_shared<GlobalMatrixT>(comm);
		}

		void clear()
		{
			if(coupling_buffer) {
				coupling_buffer->clear();
			}

			if(mass_buffer) {
				mass_buffer->clear();
			}
		}

	};

	FETransfer::FETransfer()
	{
		buffers_ = std::make_shared<Buffers>();
	}

	FETransfer::~FETransfer()
	{

	}

	void FETransfer::init(
		const std::shared_ptr<libMesh::MeshBase> &master_mesh,
		const std::shared_ptr<libMesh::DofMap>   &master_dof_map,
		const std::shared_ptr<libMesh::MeshBase> &slave_mesh,
		const std::shared_ptr<libMesh::DofMap>   &slave_dof_map,
		const std::vector<std::pair<int, int>> &tags,
		const std::vector<std::pair<int, int>> &var_pairings
		)
	{
		using namespace moonolith;

		this->master_mesh = master_mesh;
		this->slave_mesh  = slave_mesh;

		this->master_dof_map = master_dof_map;
		this->slave_dof_map  = slave_dof_map;

		this->var_pairings = var_pairings;

		predicate = std::make_shared<MasterAndSlave>();
		
		if(tags.empty()){
			auto_tag = true;
			predicate->add(0, 1);
		} else {
			auto_tag = false;
			for(auto t : tags) {   
				predicate->add(t.first, t.second);
			}
		}

		buffers_->init(slave_mesh->comm().get());
		spaces = std::make_shared<FESpacesAdapter>(master_mesh, slave_mesh, master_dof_map, slave_dof_map, var_pairings[0].first, var_pairings[0].second);
	}

	template<class TreeT>
	bool FETransfer::fill_tree(TreeT &tree)
	{
		typedef typename TreeT::DataType Adapter;

		const int n_elements_master = master_mesh->n_active_local_elem();
		const int n_elements_slave  = slave_mesh->n_active_local_elem();
		const int n_elements 		= n_elements_master + n_elements_slave;
				
		MOONOLITH_EVENT_BEGIN("create_adapters");
		
		tree.reserve(n_elements);
				
		int offset = 0;
		
		if(auto_tag){
			int space_num = 0;
			for(auto s : spaces->meshes()) {
				if(s)
				{
					bool first = true;
					libMesh::dof_id_type local_element_id = 0;
					for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it, ++local_element_id) {
						auto elem = *it;
						Adapter a(*s, elem->id(), offset + local_element_id, space_num);
						assert(!spaces->dof_map(space_num)[local_element_id].empty());
						a.set_dof_map(&spaces->dof_map(space_num)[local_element_id].global);
						tree.insert(a);
						
					}
					
					offset += s->n_active_local_elem();
				}
				
				++space_num;
			}
		} else {
			int space_num = 0;
			for(auto s : spaces->meshes()) {
				if(s) {
					bool first = true;
					
					libMesh::dof_id_type local_element_id = 0;
					for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it, ++local_element_id) {
						auto elem=*it;
						if (predicate->select(elem->subdomain_id())){
							Adapter a(*s, elem->id(), offset+local_element_id,elem->subdomain_id());
							assert(!spaces->dof_map(space_num)[local_element_id].empty());
							a.set_dof_map(&spaces->dof_map(space_num)[local_element_id].global);
							tree.insert(a);
						}
					}
					
					offset += s->n_active_local_elem();
				}
				++space_num;
			}
		}
		
		tree.root()->bound().static_bound().enlarge(1e-8);
		
		MOONOLITH_EVENT_END("create_adapters");
		return true;
	}

	template<int Dimensions>
	bool FETransfer::assemble_aux()
	{
		using namespace moonolith;
		typedef utopia::VTree<Dimensions> NTreeT;
		typedef typename NTreeT::DataContainer DataContainer;
		typedef typename NTreeT::DataType Adapter;

	
		auto tree = NTreeT::New(predicate, search_settings.max_elements, search_settings.max_depth);
		if(!fill_tree(*tree)) {
			return false;
		}
		
		return false;
	}

	bool FETransfer::assemble()
	{
		switch(slave_mesh->mesh_dimension()) {
			case 2: return assemble_aux<2>();
			case 3: return assemble_aux<3>();
			default: {
				return false;
			}
		}
	}
}
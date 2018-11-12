#include "utopia_TransferAssemblerR.hpp"

#include "utopia_LocalAssembler.hpp"
#include "utopia_BidirectionalL2LocalAssembler.hpp"

#include "utopia_Local2Global.hpp"
#include "utopia_QMortarBuilder.hpp"

#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/elem.h"

#include "utopia_libmesh.hpp"
#include "utopia_VTree.hpp"

#include "utopia_ElementDofMap.hpp"
#include "utopia_FESpacesRAdapter.hpp"

#include "moonolith_profiler.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_tree.hpp"
#include "moonolith_n_tree_mutator_factory.hpp"
#include "moonolith_n_tree_with_span_mutator_factory.hpp"
#include "moonolith_n_tree_with_tags_mutator_factory.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "par_moonolith.hpp"

#include "utopia_Socket.hpp"

#include <cmath>
#include <queue>
#include <algorithm>
#include <sstream>
#include <numeric>


/**
 * TODO:
   - Allow element-node-dof <-> node-dof convertions
   - Construction of accumulation operators
   - Discrimination of element-matrices and deletion of bad entries (e.g., duplicate volume-surface maps)
 */

namespace utopia {

	template<int Dimensions>
	class FESpaceSerializerDeserializer {
	public:
		using InputStream  = moonolith::InputStream;
		using OutputStream = moonolith::OutputStream;

		using NTreeT 		= utopia::VTree<Dimensions>;
		using DataContainer = typename NTreeT::DataContainer;
		using Adapter       = typename NTreeT::DataType;

		FESpaceSerializerDeserializer(
			const libMesh::Parallel::Communicator &comm,
			const TransferOptionsR &opts,
			const std::shared_ptr<FESpacesRAdapter> &local_spaces)
		: comm(comm),
		  m_comm(comm.get()),
		  opts(opts),
		  local_spaces(local_spaces)
		{}

		const libMesh::Parallel::Communicator &comm;
		moonolith::Communicator m_comm;
		const TransferOptionsR &opts;

		std::shared_ptr<FESpacesRAdapter> local_spaces;
		std::map<long, std::shared_ptr<FESpacesRAdapter> > spaces;
		std::map<long, std::vector<std::shared_ptr<FESpacesRAdapter> > > migrated_spaces;

		void read(
			const long ownerrank,
			const long senderrank,
			bool is_forwarding, DataContainer &data,
			InputStream &in
		) {

			CHECK_STREAM_READ_BEGIN("vol_proj", in);

			std::shared_ptr<FESpacesRAdapter> proc_space = std::make_shared<FESpacesRAdapter>(m_comm);

			read_spaces(in, *proc_space, comm, comm);

			if (!is_forwarding) {
				assert(!spaces[ownerrank]);
				spaces[ownerrank] = proc_space;
			} else {
				migrated_spaces[ownerrank].push_back(proc_space);
			}

			data.reserve(data.size() + 3000);

			long offset = 0;

			if(opts.tags.empty()){
				int space_num = 0;
				for(auto s : proc_space->spaces()) {
					if(s) {
						//ID_FIX this should be fine n_elem is actually local sence the mesh is a SerialMesh
						for (int i=0; i<s->n_elem(); i++) {
							data.push_back(Adapter(*s, i, offset + i,space_num));
							assert(!proc_space->dof_map(space_num)[i].empty());
							assert(!proc_space->dof_map_reverse(space_num)[i].empty());
							data.back().set_dof_map(&proc_space->dof_map(space_num)[i].global);
							data.back().set_dof_map_reverse(&proc_space->dof_map_reverse(space_num)[i].global);
						}

						offset += s->n_elem();

					}

					++space_num;
				}
			} else {
				int space_num = 0;
				for(auto s : proc_space->spaces()) {
					if(s) {
						for (int i=0; i<s->n_elem(); i++) {
							const libMesh::Elem * elem = s->elem_ptr(i);
							//Volume Tag
							int volume_tag = elem->subdomain_id();
							data.push_back(Adapter(*s, i, offset + i, volume_tag));
							assert(!proc_space->dof_map(space_num)[i].empty());
							assert(!proc_space->dof_map_reverse(space_num)[i].empty());
							data.back().set_dof_map(&proc_space->dof_map(space_num)[i].global);
							data.back().set_dof_map_reverse(&proc_space->dof_map_reverse(space_num)[i].global);
						}

						offset += s->n_elem();

					}

					++space_num;
				}
			}

			CHECK_STREAM_READ_END("vol_proj", in);
		};

		void write(
			const long ownerrank, const long recvrank,
			const std::vector<long>::const_iterator &begin,
			const std::vector<long>::const_iterator &end,
			const DataContainer &data,
			OutputStream &out) {

			CHECK_STREAM_WRITE_BEGIN("vol_proj", out);

			if (ownerrank == m_comm.rank()) {
				write_element_selection(begin, end, *local_spaces, out);
			} else {
				auto it = spaces.find(ownerrank);
				assert(it != spaces.end());
				std::shared_ptr<FESpacesRAdapter> spaceptr = it->second;
				assert(std::distance(begin, end) > 0);
				write_element_selection(begin, end, *spaceptr, out);
			}

			CHECK_STREAM_WRITE_END("vol_proj", out);
		};

	};

	template class FESpaceSerializerDeserializer<2>;
	template class FESpaceSerializerDeserializer<3>;

	template<int Dimensions>
	class DefaultAlgorithmR final : public TransferAssemblerR::AlgorithmR {
	public:
		using FunctionSpace = utopia::LibMeshFunctionSpace;
		using SparseMatrix  = utopia::USparseMatrix;
		using MeshBase      = libMesh::MeshBase;
		using DofMap        = libMesh::DofMap;
		using NTreeT 		= utopia::VTree<Dimensions>;
		using DataContainer = typename NTreeT::DataContainer;
		using Adapter       = typename NTreeT::DataType;

		DefaultAlgorithmR(
			const std::shared_ptr<MeshBase> &from_mesh,
			const std::shared_ptr<DofMap>   &from_dofs,
			const std::shared_ptr<DofMap>   &from_dofs_r,
			const std::shared_ptr<MeshBase> &to_mesh,
			const std::shared_ptr<DofMap>   &to_dofs,
			const std::shared_ptr<DofMap>   &to_dofs_r,
			const TransferOptionsR &opts,
			const std::shared_ptr<LocalAssembler> &assembler,
			const std::shared_ptr<Local2Global>  &local2global)
		{
			init(from_mesh, from_dofs, from_dofs_r, to_mesh, to_dofs, to_dofs_r, opts, assembler, local2global);
		}

		void init(
			const std::shared_ptr<MeshBase> &from_mesh,
			const std::shared_ptr<DofMap>   &from_dofs,
			const std::shared_ptr<DofMap>   &from_dofs_r,
			const std::shared_ptr<MeshBase> &to_mesh,
			const std::shared_ptr<DofMap>   &to_dofs,
			const std::shared_ptr<DofMap>   &to_dofs_r,
			const TransferOptionsR &opts,
			const std::shared_ptr<LocalAssembler> &assembler,
			const std::shared_ptr<Local2Global>   &local2global)
		{
			this->from_mesh = from_mesh;
			this->from_dofs = from_dofs;
			this->from_dofs_r = from_dofs_r;
			this->to_mesh 	= to_mesh;
			this->to_dofs	= to_dofs;
			this->to_dofs_r	= to_dofs_r;
			this->opts 		= opts;

			this->assembler = assembler;
			this->local2global = local2global;

			this->comm = moonolith::Communicator(from_mesh->comm().get());
			std::cout<<"FESpacesRAdapte"<<std::endl;
			this->local_spaces = std::make_shared<FESpacesRAdapter>(from_mesh, to_mesh, from_dofs, to_dofs, from_dofs_r, to_dofs_r, opts.from_var_num, opts.to_var_num, opts.from_var_num_r, opts.to_var_num_r);
		   
			predicate = std::make_shared<moonolith::MasterAndSlave>();
            std::cout<<"FESpacesRAdapte after"<<std::endl;
			if(opts.tags.empty()){
				predicate->add(0, 1);
			} else {
				for(auto t : opts.tags) {
					predicate->add(t.first, t.second);
				}
			}
			std::cout<<"FESpacesRAdapte after 2"<<std::endl;
		}

		void pre_assemble()
		{
			const std::size_t n_forms = mat_buffer.size();

			std::cout<<"pre_assemble"<<std::endl;
			
			for(std::size_t i = 0; i < n_forms; ++i) {
				mat_buffer[i] = std::make_shared< moonolith::SparseMatrix<double> >(comm);
				local_element_matrices_sum[i] = 0.;

				switch(assembler->type(i)) {
					
					case LocalAssembler::MASTER_X_SLAVE: 
					{
						
					    std::cout<<"MASTER_X_SLAVE"<<std::endl;
						mat_buffer[i]->set_size(to_dofs->n_dofs(), from_dofs->n_dofs());
						break;
					}

					case LocalAssembler::SLAVE_X_SLAVE:
					{
						std::cout<<"SLAVE_X_SLAVE"<<std::endl;
						mat_buffer[i]->set_size(to_dofs->n_dofs(), to_dofs->n_dofs());
						break;
					}

					case LocalAssembler::MASTER_X_MASTER:
					{
						std::cout<<"MASTER_X_MASTER"<<std::endl;

						mat_buffer[i]->set_size(from_dofs_r->n_dofs(), from_dofs_r->n_dofs());
						break;
					}

					case LocalAssembler::SLAVE_X_MASTER: 
					{
						
                        std::cout<<"SLAVE_X_MASTER"<<std::endl;
						mat_buffer[i]->set_size(from_dofs_r->n_dofs(), to_dofs_r->n_dofs());
						break;
					}

					default:
					{
						assert(false);
						break;
					}
				}
			}
		}

		bool assemble(Adapter &master,
					  Adapter &slave)

		{
			

		    std::cout<<"assemble"<<std::endl;

			//FIXME assuming elements are all the same
		 	auto master_type = from_dofs->variable(opts.from_var_num).type();
		 	auto slave_type  = to_dofs->variable(opts.to_var_num).type();

			const auto &master_mesh = master.space();;
			const auto &slave_mesh  = slave.space();

			const int src_index  = master.element();
			const int dest_index = slave.element();

			auto &master_el = *master_mesh.elem(src_index);
			auto &slave_el  = *slave_mesh.elem(dest_index);

			for(auto &mat_i : elemmat) {
				mat_i.zero();
			}

			if(assembler->assemble(master_el, master_type, slave_el, slave_type, elemmat)) {
				
				for(std::size_t i = 0; i < elemmat.size(); ++i) {	
					auto &mat_i = elemmat[i];
					auto partial_sum = std::accumulate(mat_i.get_values().begin(), mat_i.get_values().end(), libMesh::Real(0.0));
					local_element_matrices_sum[i] += partial_sum;

					switch(assembler->type(i)) {
						
						case LocalAssembler::MASTER_X_SLAVE: 
						{
							const auto &master_dofs = master.dof_map();
							const auto &slave_dofs  = slave.dof_map();

							local2global->apply(master_dofs, slave_dofs, elemmat[i], *mat_buffer[i]);
							break;
						}

						case LocalAssembler::SLAVE_X_SLAVE:
						{
							const auto &slave_dofs  = slave.dof_map();

							local2global->apply(slave_dofs, slave_dofs, elemmat[i], *mat_buffer[i]);
							break;
						}

						case LocalAssembler::MASTER_X_MASTER:
						{
							
							const auto &master_dofs_r = master.dof_map_reverse();

							local2global->apply(master_dofs_r, master_dofs_r, elemmat[i], *mat_buffer[i]);
							break;
						}

						case LocalAssembler::SLAVE_X_MASTER: 
						{
							
							const auto &master_dofs_r = master.dof_map_reverse();
							const auto &slave_dofs_r  = slave.dof_map_reverse();

							local2global->apply(slave_dofs_r, master_dofs_r, elemmat[i], *mat_buffer[i]);
							break;
						}

						default:
						{
							assert(false);
							break;
						}
					}
				}

				return true;
			} else {
				return false;
			}
		}

		void print_stats()
		{
			double total_intersection_volume = 0.;
			{
				auto l2_assembler = std::dynamic_pointer_cast<BidirectionalL2LocalAssembler>(assembler);
				if(l2_assembler) {
					total_intersection_volume = l2_assembler->get_q_builder().get_total_intersection_volume();

					double volumes[2] = { local_element_matrices_sum[0], total_intersection_volume };
					comm.all_reduce(volumes, 2, moonolith::MPISum());

					if(comm.is_root()) {
						std::cout << "sum(B): " << volumes[0] << ", vol(I): " << volumes[1] << std::endl;
					}
				}
			}
		}

		void post_assemble(std::size_t buffer_num)
		{
			SparseMatrix &mat = *mats_[buffer_num];

			libMesh::dof_id_type n_dofs_on_proc_trial = 0;
			libMesh::dof_id_type n_dofs_on_proc_test  = 0;

			switch(assembler->type(buffer_num)) {
				
				case LocalAssembler::MASTER_X_SLAVE: 
				{
					n_dofs_on_proc_trial = from_dofs->n_local_dofs();
					n_dofs_on_proc_test  = to_dofs->n_local_dofs();
					break;
				}

				case LocalAssembler::SLAVE_X_SLAVE:
				{
					n_dofs_on_proc_trial = to_dofs->n_local_dofs();
					n_dofs_on_proc_test  = to_dofs->n_local_dofs();
					break;
				}

				case LocalAssembler::MASTER_X_MASTER:
				{
					n_dofs_on_proc_trial = from_dofs_r->n_local_dofs();
					n_dofs_on_proc_test  = from_dofs_r->n_local_dofs();
					break;
				}

				case LocalAssembler::SLAVE_X_MASTER: 
				{
					n_dofs_on_proc_trial = to_dofs_r->n_local_dofs();
					n_dofs_on_proc_test  = from_dofs_r->n_local_dofs();
					break;
				}

				default:
				{
					assert(false);
					break;
				}
			}

			local2global->redistribute(comm, n_dofs_on_proc_trial, n_dofs_on_proc_test, *mat_buffer[buffer_num]);

			SizeType m_max_row_entries = mat_buffer[buffer_num]->local_max_entries_x_col();
			comm.all_reduce(&m_max_row_entries, 1, moonolith::MPIMax());

			USparseMatrix mat_x = utopia::local_sparse(n_dofs_on_proc_test, n_dofs_on_proc_trial, m_max_row_entries);

			{
				utopia::Write<utopia::USparseMatrix> write(mat_x);
				for (auto it = mat_buffer[buffer_num]->iter(); it; ++it) {
					mat_x.set(it.row(), it.col(), *it);

				}
			}

			if(opts.n_var == 1) {
				mat = std::move(mat_x);
				return;
			}

			auto s_mat_x = local_size(mat_x);
			mat = local_sparse(s_mat_x.get(0), s_mat_x.get(1), opts.n_var * m_max_row_entries);

			utopia::Write<USparseMatrix> w_mat(mat);
			utopia::each_read(mat_x, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
				for(utopia::SizeType d = 0; d < opts.n_var; ++d) {
					mat.set(i + d, j + d, value);
				}
			});
		}

		bool assemble(std::vector<std::shared_ptr<SparseMatrix> > &mats)
		{
			std::cout<<" I am  in assemble" <<std::endl;

			if(assembler->n_forms() != mats.size()) {
				mats.resize(assembler->n_forms());
				std::cout<<"assembler->n_forms()"<<assembler->n_forms()<<std::endl;
			}

			for(auto &mat_ptr : mats) {
				if(!mat_ptr) {
					mat_ptr = std::make_shared<SparseMatrix>();
				}
			}

			init_buffers(assembler->n_forms());
			mats_ = mats;
			return assemble_aux();
		}

		bool assemble_aux()
		{
			using namespace moonolith;

			init_tree();

			std::map<long, std::shared_ptr<FESpacesRAdapter> > spaces;
			std::map<long, std::vector<std::shared_ptr<FESpacesRAdapter> > > migrated_spaces;

			FESpaceSerializerDeserializer<Dimensions> serializer(
				from_mesh->comm(),
				opts,
				local_spaces);

			auto read = [&serializer] (
				const long ownerrank,
				const long senderrank,
				bool is_forwarding, DataContainer &data,
				InputStream &in
			) {
				serializer.read(ownerrank, senderrank, is_forwarding, data, in);
			};

			auto write = [&serializer] (
				const long ownerrank, const long recvrank,
				const std::vector<long>::const_iterator &begin,
				const std::vector<long>::const_iterator &end,
				const DataContainer &data,
				OutputStream &out
			) {
				serializer.write(ownerrank, recvrank, begin, end, data, out);
			};

			long n_false_positives = 0, n_intersections = 0;

			auto fun = [&](Adapter &master, Adapter &slave) -> bool {
				if(this->assemble(master, slave)) {
					n_intersections++;
					return true;
				} else {
					n_false_positives++;
					return false;
				}
			};

			pre_assemble();

			moonolith::search_and_compute(comm, tree, predicate, read, write, fun, settings);

			print_stats();
			
			for(std::size_t i = 0; i < mats_.size(); ++i) {
				post_assemble(i);
			}

			long n_total_candidates = n_intersections + n_false_positives;
			long n_collection[3] = {n_intersections, n_total_candidates, n_false_positives};

			comm.all_reduce(n_collection, 3, moonolith::MPISum());

			if (comm.is_root()) {
				std::cout << "n_intersections: " << n_collection[0]
				<< ", n_total_candidates: " 	 << n_collection[1]
				<< ", n_false_positives: " 	     << n_collection[2] << std::endl;
			}

			return true;
		}

		void init_tree()
		{
			using namespace moonolith;

			const auto n_elements_from = from_mesh->n_active_local_elem();
			const auto n_elements_to   = to_mesh->n_active_local_elem();
			const auto n_elements 	  = n_elements_from + n_elements_to;

			MOONOLITH_EVENT_BEGIN("create_adapters");

			tree = NTreeT::New(predicate, settings.max_elements, settings.max_depth);
			tree->reserve(n_elements);

			int offset = 0;
			if(opts.tags.empty()){
				int space_num = 0;

				for(auto s : local_spaces->spaces()) {

					if(s)
					{
						bool first = true;
						libMesh::dof_id_type local_element_id = 0;
						for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it, ++local_element_id) {
							auto elem = *it;
							Adapter a(*s, elem->id(), offset+local_element_id,space_num);
							assert(!local_spaces->dof_map(space_num)[local_element_id].empty());
							assert(!local_spaces->dof_map_reverse(space_num)[local_element_id].empty());
							a.set_dof_map(&local_spaces->dof_map(space_num)[local_element_id].global);
							a.set_dof_map_reverse(&local_spaces->dof_map_reverse(space_num)[local_element_id].global);
							tree->insert(a);

						}

						offset += s->n_active_local_elem();
					}

					++space_num;
				}

			} else {

				int space_num = 0;
				for(auto s : local_spaces->spaces()) {
					if(s) {

						bool first = true;

						libMesh::dof_id_type local_element_id = 0;
						for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it, ++local_element_id) {
							auto elem=*it;
							if (predicate->select(elem->subdomain_id())){
								Adapter a(*s, elem->id(), offset+local_element_id,elem->subdomain_id());
								assert(!local_spaces->dof_map(space_num)[local_element_id].empty());
								assert(!local_spaces->dof_map_reverse(space_num)[local_element_id].empty());
								a.set_dof_map(&local_spaces->dof_map(space_num)[local_element_id].global);
								a.set_dof_map_reverse(&local_spaces->dof_map_reverse(space_num)[local_element_id].global);
								tree->insert(a);
							}
						}

						offset += s->n_active_local_elem();
					}

					++space_num;
				}
			}

			tree->root()->bound().static_bound().enlarge(1e-8);

			MOONOLITH_EVENT_END("create_adapters");
		}

		void init_buffers(const SizeType n)
		{
			std::cout<<"init_buffers"<<n<<std::endl;
			mat_buffer.resize(n);
			elemmat.resize(n);
			local_element_matrices_sum.resize(n);
		}

	private:
		std::shared_ptr<MeshBase> from_mesh;
		std::shared_ptr<DofMap>   from_dofs;
		std::shared_ptr<DofMap>   from_dofs_r;
		std::shared_ptr<MeshBase> to_mesh;
		std::shared_ptr<DofMap>   to_dofs;
		std::shared_ptr<DofMap>   to_dofs_r;
		TransferOptionsR opts;

		std::shared_ptr<LocalAssembler> assembler;
		std::shared_ptr<Local2Global> local2global;

		moonolith::Communicator comm;
		moonolith::SearchSettings settings;
		std::shared_ptr<moonolith::MasterAndSlave> predicate;

		std::shared_ptr<NTreeT> tree;
		std::shared_ptr<FESpacesRAdapter> local_spaces;

		std::vector< libMesh::DenseMatrix<libMesh::Real> > elemmat;
		std::vector< libMesh::Real > local_element_matrices_sum;

		std::vector< std::shared_ptr< moonolith::SparseMatrix<double> > > mat_buffer;
		std::vector<std::shared_ptr<SparseMatrix>> mats_;
	};

	template class DefaultAlgorithmR<2>;
	template class DefaultAlgorithmR<3>;


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	TransferAssemblerR::TransferAssemblerR(
		const std::shared_ptr<LocalAssembler> &assembler,
		const std::shared_ptr<Local2Global> &local2global)
	: assembler_(assembler), local2global_(local2global)
	{}

	TransferAssemblerR::~TransferAssemblerR() {}


	bool TransferAssemblerR::assemble(
		const std::shared_ptr<MeshBase> &from_mesh,
		const std::shared_ptr<DofMap>   &from_dofs,
		const std::shared_ptr<DofMap>   &from_dofs_r,
		const std::shared_ptr<MeshBase> &to_mesh,
		const std::shared_ptr<DofMap>   &to_dofs,
		const std::shared_ptr<DofMap>   &to_dofs_r,
		std::vector<std::shared_ptr<SparseMatrix> > &mats,
		const TransferOptionsR &opts)
	{
		assert(assembler_    && "assembler is required");
		assert(local2global_ && "local2global");

		///////////////////////////

		moonolith::Communicator comm(from_mesh->comm().get());

		if(Utopia::instance().verbose()) {
			moonolith::root_describe("---------------------------------------\n"
				"begin: utopia::TransferAssembler::assemble",
				comm, std::cout);
		}

		Chrono c;
		c.start();

		///////////////////////////

		if(from_mesh->mesh_dimension() == 2) {
			std::cout<<"I am here"<<std::endl;
			algorithm_ = std::make_shared<DefaultAlgorithmR<2>>(from_mesh, from_dofs, from_dofs_r, to_mesh, to_dofs, to_dofs_r, opts, assembler_, local2global_);
		} else if(from_mesh->mesh_dimension() == 3) {
			algorithm_ = std::make_shared<DefaultAlgorithmR<3>>(from_mesh, from_dofs, from_dofs_r, to_mesh, to_dofs, to_dofs_r, opts, assembler_, local2global_);
		} else {
			assert(false && "dimension not supported");
			return false;
		}

		std::cout<<"I am in"<<std::endl;

		bool ok = algorithm_->assemble(mats);

		std::cout<<"I am out"<<std::endl;


		///////////////////////////
		c.stop();

		if(Utopia::instance().verbose()) {
			std::stringstream ss;
			ss << "end: utopia::TransferAssembler::assemble\n";
			ss << c;
			ss << "---------------------------------------";
			moonolith::root_describe(ss.str(), comm, std::cout);
		}
		///////////////////////////

		return ok;
	}

	bool TransferAssemblerR::assemble(
            const std::shared_ptr<MeshBase> &from_mesh,
			const std::shared_ptr<DofMap>   &from_dofs,
			const std::shared_ptr<DofMap>   &from_dofs_r,
			const std::shared_ptr<MeshBase> &to_mesh,
			const std::shared_ptr<DofMap>   &to_dofs,
			const std::shared_ptr<DofMap>   &to_dofs_r,
			SparseMatrix &B,
			const TransferOptionsR &opts
		)
	{
		std::vector<std::shared_ptr<SparseMatrix> > mats;
		mats.push_back(make_ref(B));
		std::cout<<"I am in TransferAssemblerR::assemble"<<std::endl;
		return assemble(from_mesh, from_dofs, from_dofs_r, to_mesh, to_dofs, to_dofs_r, mats, opts);
	}
}

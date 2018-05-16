#include "utopia_assemble_volume_transfer.hpp"

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
#include "par_moonolith.hpp"

#include "utopia_Socket.hpp"

#include <cmath>
#include <queue>
#include <algorithm>
#include <sstream>
#include <numeric>


namespace utopia {

	static void assemble_biorth_weights(
		const libMesh::Elem &el,
		const int dim,
		const libMesh::FEType &var_type,
		const int el_order,
		libMesh::DenseMatrix<libMesh::Real> &weights)
	{
		std::unique_ptr<libMesh::FEBase> biorth_elem = libMesh::FEBase::build(dim, var_type);

		const int order = order_for_l2_integral(dim, el, el_order, el, el_order);

		libMesh::QGauss qg(dim, libMesh::Order(order));
		biorth_elem->attach_quadrature_rule(&qg);
		biorth_elem->reinit(&el);
		mortar_assemble_weights(*biorth_elem, weights);
	}

	class Local2Global {
	public:
		using LocalMatrix = libMesh::DenseMatrix<libMesh::Real>;

		Local2Global(const bool use_set_instead_of_add)
		: use_set_instead_of_add(use_set_instead_of_add)
		{}

		void apply(
			const std::vector<long> &trial,
			const std::vector<long> &test,
			const LocalMatrix &local_mat,
			moonolith::SparseMatrix<double> &global_mat) const
		{
			if(use_set_instead_of_add) {
				set(trial, test, local_mat, global_mat);
			} else {
				add(trial, test, local_mat, global_mat);
			}
		}

		void add(
			const std::vector<long> &trial,
			const std::vector<long> &test,
			const LocalMatrix &local_mat,
			moonolith::SparseMatrix<double> &global_mat) const
		{
			for(std::size_t i = 0; i < test.size(); ++i) {
				const auto dof_I = test[i];

				for(std::size_t j = 0; j < trial.size(); ++j) {
					const auto dof_J = trial[j];

					global_mat.add(dof_I, dof_J, local_mat(i, j));
				}
			}
		}

		void set(
			const std::vector<long> &trial,
			const std::vector<long> &test,
			const LocalMatrix &local_mat,
			moonolith::SparseMatrix<double> &global_mat) const
		{
			for(std::size_t i = 0; i < test.size(); ++i) {
				const auto dof_I = test[i];

				for(std::size_t j = 0; j < trial.size(); ++j) {
					const auto dof_J = trial[j];

					global_mat.set(dof_I, dof_J, local_mat(i, j));
				}
			}
		}

		void redistribute(
			moonolith::Communicator &comm,
			const long n_local_dofs_trial,
			const long n_local_dofs_test,
			moonolith::SparseMatrix<double> &global_mat)
		{
			std::vector<moonolith::Integer> range_master(comm.size() + 1, 0);
			std::vector<moonolith::Integer> range_slave(comm.size()  + 1, 0);

			range_master[comm.rank() + 1] += static_cast<unsigned int>(n_local_dofs_trial);
			range_slave [comm.rank() + 1] += static_cast<unsigned int>(n_local_dofs_test);

			comm.all_reduce(&range_master[0], range_master.size(), moonolith::MPISum());
			comm.all_reduce(&range_slave[0],  range_slave.size(),  moonolith::MPISum());

			std::partial_sum(range_master.begin(), range_master.end(), range_master.begin());
			std::partial_sum(range_slave.begin(),  range_slave.end(),  range_slave.begin());

			moonolith::Redistribute< moonolith::SparseMatrix<double> > redist(comm.get_mpi_comm());

			if(use_set_instead_of_add) {
				redist.apply(range_slave, global_mat, moonolith::Assign<double>());
			} else {
				redist.apply(range_slave, global_mat, moonolith::AddAssign<double>());
			}

			assert(range_slave.empty() == range_master.empty() || range_master.empty());
		}

		bool use_set_instead_of_add;
	};

	class LocalAssembler {
	public:
		using Elem = libMesh::Elem;
		using FEType = libMesh::FEType;
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;

		virtual ~LocalAssembler() {}
		virtual bool assemble(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			Matrix &mat
			) = 0;
	};

	class QMortarBuilder {
	public:
		using Elem = libMesh::Elem;
		using FEType = libMesh::FEType;

		virtual ~QMortarBuilder() {}

		virtual bool build(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			QMortar &q_trial,
			QMortar &q_test) = 0;

		virtual double get_total_intersection_volume() const = 0;

	};

	class QMortarBuilder2 final : public QMortarBuilder {
	public:
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;
		using Real = libMesh::Real;

		QMortarBuilder2()
		: total_intersection_volume(0), composite_ir(2)
		{}

		virtual bool build(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			QMortar &q_trial,
			QMortar &q_test) override
		{
			make_polygon(trial, trial_pts);
			make_polygon(test, test_pts);

			if(intersect_2D(trial_pts, test_pts, intersection)) {
				total_intersection_volume += fabs(isector.polygon_area_2(intersection.m(), &intersection.get_values()[0]));
				const libMesh::Real weight = isector.polygon_area_2(test_pts.m(), &test_pts.get_values()[0]);

				const int order = order_for_l2_integral(2, trial, trial_type.order, test, test_type.order);
				make_composite_quadrature_2D(intersection, weight, order, composite_ir);
				auto trial_trans  = std::make_shared<AffineTransform2>(trial);
				auto test_trans   = std::make_shared<AffineTransform2>(test);

				transform_to_reference(*trial_trans,  trial.type(), composite_ir, q_trial);


				// if(vol2surf) {
				// 	transform_to_reference_surf(*dest_trans, dest_el.type(), composite_ir, dest_ir);
				// } else {
				transform_to_reference(*test_trans,   test.type(),  composite_ir, q_test);
				// }


				return true;
			} else {
				return false;
			}
		}

		double get_total_intersection_volume() const override
		{
			return total_intersection_volume;
		}

	private:
		Matrix trial_pts;
		Matrix test_pts;
		Matrix intersection;

		Real total_intersection_volume;
		QMortar composite_ir;
		Intersector isector;
	};


	class QMortarBuilder3 final : public QMortarBuilder {
	public:
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;
		using Real = libMesh::Real;

		QMortarBuilder3()
		: total_intersection_volume(0), composite_ir(3)
		{}

		virtual bool build(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			QMortar &q_trial,
			QMortar &q_test) override
		{
			make_polyhedron(trial, trial_poly);
			make_polyhedron(test,  test_poly);

			if(intersect_3D(trial_poly, test_poly, intersection)) {
				total_intersection_volume += compute_volume(intersection);
				const libMesh::Real weight = compute_volume(test_poly);
				auto trial_trans = std::make_shared<AffineTransform3>(trial);
				std::shared_ptr<Transform> test_trans;

				bool vol2surf = false;

				if(is_tri(test.type()) || is_quad(test.type())) {
					test_trans = std::make_shared<Transform2>(test);
					vol2surf = true;
				} else {
					test_trans = std::make_shared<AffineTransform3>(test);
				}

				const int order = order_for_l2_integral(2, trial, trial_type.order, test, test_type.order);

				if(vol2surf) {
					shell_poly.resize(intersection.n_nodes, 3);
					std::copy(intersection.points, intersection.points + intersection.n_nodes * intersection.n_dims, &shell_poly.get_values()[0]);
					make_composite_quadrature_on_surf_3D(shell_poly, 1./weight, order, composite_ir);

				// plot_polygon(3, test_poly.n_nodes, test_poly.points, "polygon/" + std::to_string(comm.rank()) + "/p");
				} else {
					make_composite_quadrature_3D(intersection, weight, order, composite_ir);
				}

				transform_to_reference(*trial_trans, trial.type(), composite_ir,  q_trial);

				if(vol2surf) {
					transform_to_reference_surf(*test_trans, test.type(), composite_ir, q_test);
				} else {
					transform_to_reference(*test_trans, test.type(), composite_ir,  q_test);
				}

				return true;
			} else {
				return false;
			}
		}

		double get_total_intersection_volume() const override
		{
			return total_intersection_volume;
		}

	private:
		Real total_intersection_volume;
		QMortar composite_ir;
		Intersector isector;

		Polyhedron trial_poly, test_poly;
		Polyhedron intersection, temp_poly;

		Matrix shell_poly;
	};

	class L2LocalAssembler final : public LocalAssembler {
	public:
		L2LocalAssembler(const int dim, const bool use_biorth)
		: dim(dim),
		  use_biorth(use_biorth),
		  must_compute_biorth(use_biorth),
		  composite_ir(dim),
		  q_trial(dim),
		  q_test(dim)
		{
			if(dim == 2) {
				q_builder = std::make_shared<QMortarBuilder2>();
			} else {
				assert(dim == 3);
				q_builder = std::make_shared<QMortarBuilder3>(); 
			}
		}

		virtual bool assemble(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			Matrix &mat
			) override
		{
			auto trial_fe   = libMesh::FEBase::build(trial.dim(), trial_type);
			auto test_fe    = libMesh::FEBase::build(test.dim(),  test_type);
			const int order = order_for_l2_integral(dim, trial, trial_type.order, test, test_type.order);

			if(!q_builder->build(trial, trial_type, test, test_type, q_trial, q_test)) {
				return false;
			}

			init_biorth(test, test_type);
			init_fe(trial, trial_type, test, test_type);

			trial_fe->attach_quadrature_rule(&q_trial);
			trial_fe->get_phi();
			trial_fe->reinit(&trial);

			test_fe->attach_quadrature_rule(&q_test);
			test_fe->get_phi();
			test_fe->get_JxW();
			test_fe->reinit(&test);

			if(use_biorth) {
				mortar_assemble_weighted_biorth(*trial_fe, *test_fe, biorth_weights, mat);
			} else {
				mortar_assemble(*trial_fe, *test_fe, mat);
			}

			return true;
		}

		void init_fe(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type)
		{
			if(trial_fe) return;

			trial_fe = libMesh::FEBase::build(trial.dim(), trial_type);
			test_fe  = libMesh::FEBase::build(test.dim(),  test_type);
		}

		void init_biorth(const Elem &test, FEType test_type)
		{
			if(!use_biorth) return;

			assemble_biorth_weights(
				test,
				test.dim(),
				test_type,
				test_type.order,
				biorth_weights);

			must_compute_biorth = false;
		}


		inline const QMortarBuilder &get_q_builder() const
		{
			assert(q_builder);
			return *q_builder;
		}
	private:
		int dim;
		bool use_biorth;
		bool must_compute_biorth;
		QMortar composite_ir;
		QMortar q_trial;
		QMortar q_test;

		Matrix biorth_weights;

		std::shared_ptr<QMortarBuilder> q_builder;
		std::unique_ptr<libMesh::FEBase> trial_fe, test_fe;
	};

	class InterpolationLocalAssembler final : public LocalAssembler {
	public:
		virtual bool assemble(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			Matrix &mat
			) override
		{
			return false;
		}
	};
}

using namespace libMesh;

namespace utopia {
	
	
	template<int Dimensions, class Fun>
	static bool Assemble(
		const std::shared_ptr<MeshBase> &master,
		const std::shared_ptr<MeshBase> &slave,
		const std::shared_ptr<DofMap> &dof_master,
		const std::shared_ptr<DofMap> &dof_slave,
		const unsigned int &from_var_num,
		const unsigned int &to_var_num,
		Fun process_fun,
		const moonolith::SearchSettings &settings,
		bool use_biorth,
		const std::vector< std::pair<int, int> > &tags,
		int n_var)
	{
		
		
		using namespace moonolith;
		
		typedef VTree<Dimensions> NTreeT;
		typedef typename NTreeT::DataContainer DataContainer;
		typedef typename NTreeT::DataType Adapter;
		
		const long maxNElements = settings.max_elements;
		const long maxDepth = settings.max_depth;
		
		
		const auto &master_mesh = master;
		const auto &slave_mesh  = slave;
		const int n_elements_master = master_mesh->n_active_local_elem();
		const int n_elements_slave  = slave_mesh->n_active_local_elem();
		const int n_elements 		= n_elements_master + n_elements_slave;
		
		
		const Parallel::Communicator &libmesh_comm_master = master_mesh->comm();
		const Parallel::Communicator &libmesh_comm_slave = slave_mesh->comm();

		moonolith::Communicator comm(libmesh_comm_master.get());
		
		
		auto predicate = std::make_shared<MasterAndSlave>();
		
		if(tags.empty()){
			//           std::cout<<"tags.empty()==>"<<tags.empty()<<std::endl;
			predicate->add(0, 1);
		}
		else{
			
			for(auto t : tags)
				{   predicate->add(t.first, t.second);
				//                std::cout<<"t.first==>"<<t.first<<std::endl;
				//                std::cout<<"t.second==>"<<t.second<<std::endl;
				}
			}


			MOONOLITH_EVENT_BEGIN("create_adapters");

			std::shared_ptr<NTreeT> tree = NTreeT::New(predicate, maxNElements, maxDepth);
			tree->reserve(n_elements);


			auto local_spaces = std::make_shared<FESpacesAdapter>(master, slave, dof_master, dof_slave, from_var_num, to_var_num);

			int offset = 0;


			if(tags.empty()){

				int space_num = 0;
				for(auto s : local_spaces->spaces()) {
					if(s)
					{
						bool first = true;
						libMesh::dof_id_type local_element_id = 0;
						for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it, ++local_element_id) {
							auto elem=*it;
							Adapter a(*s, elem->id(), offset+local_element_id,space_num);
							assert(!local_spaces->dof_map(space_num)[local_element_id].empty());
							a.set_dof_map(&local_spaces->dof_map(space_num)[local_element_id].global);
							tree->insert(a);

						}

						offset += s->n_active_local_elem();
					}

					++space_num;
				}


			}
			else
			{
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
								a.set_dof_map(&local_spaces->dof_map(space_num)[local_element_id].global);
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


			std::map<long, std::shared_ptr<FESpacesAdapter> > spaces;
			std::map<long, std::vector<std::shared_ptr<FESpacesAdapter> > > migrated_spaces;


			auto read = [&spaces, &migrated_spaces, comm, &libmesh_comm_master, &libmesh_comm_slave, &tags]
			(
				const long ownerrank,
				const long senderrank,
				bool is_forwarding, DataContainer &data,
				InputStream &in
				) {
				CHECK_STREAM_READ_BEGIN("vol_proj", in);

				std::shared_ptr<FESpacesAdapter> proc_space = std::make_shared<FESpacesAdapter>(comm);

				read_spaces(in, *proc_space, libmesh_comm_master, libmesh_comm_slave);

				if (!is_forwarding) {
					assert(!spaces[ownerrank]);
					spaces[ownerrank] = proc_space;
				} else {
					migrated_spaces[ownerrank].push_back(proc_space);
				}

				data.reserve(data.size() + 3000);



				long offset = 0;

				if(tags.empty()){
					int space_num = 0;
					for(auto s : proc_space->spaces()) {
						if(s) {
						//ID_FIX this should be fine n_elem is actually local sence the mesh is a SerialMesh
							for (int i=0; i<s->n_elem(); i++) {
								data.push_back(Adapter(*s, i, offset + i,space_num));
								assert(!proc_space->dof_map(space_num)[i].empty());
								data.back().set_dof_map(&proc_space->dof_map(space_num)[i].global);
							}

							offset += s->n_elem();

						}

						++space_num;
					}
				}

				else{
					int space_num = 0;
					for(auto s : proc_space->spaces()) {
						if(s) {
							for (int i=0; i<s->n_elem(); i++) {
								const Elem * elem = s->elem_ptr(i);
							//Volume Tag
								int volume_tag = elem->subdomain_id();
								data.push_back(Adapter(*s, i, offset + i,volume_tag) );
								assert(!proc_space->dof_map(space_num)[i].empty());
								data.back().set_dof_map(&proc_space->dof_map(space_num)[i].global);
							}

							offset += s->n_elem();

						}

						++space_num;
					}
				}

				CHECK_STREAM_READ_END("vol_proj", in);
			};


			auto write = [&local_spaces, &spaces, &comm] (
				const long ownerrank, const long recvrank,
				const std::vector<long>::const_iterator &begin,
				const std::vector<long>::const_iterator &end,
				const DataContainer &data,
				OutputStream &out) {

				CHECK_STREAM_WRITE_BEGIN("vol_proj", out);


				if (ownerrank == comm.rank()) {

					write_element_selection(begin, end, *local_spaces, out);


				} else {

					auto it = spaces.find(ownerrank);
					assert(it != spaces.end());
					std::shared_ptr<FESpacesAdapter> spaceptr = it->second;
					assert(std::distance(begin, end) > 0);
					write_element_selection(begin, end, *spaceptr, out);

				}


				CHECK_STREAM_WRITE_END("vol_proj", out);

			};


			long n_false_positives = 0, n_intersections = 0;


			auto fun = [&n_false_positives, &n_intersections, &process_fun](Adapter &master, Adapter &slave) -> bool {
				bool ok = process_fun(master, slave);

				if(ok) {
					n_intersections++;
					return true;
				} else {
					n_false_positives++;
					return false;
				}

				return true;
			};

			moonolith::search_and_compute(comm, tree, predicate, read, write, fun, settings);

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

		static void scale_polyhedron(const double scaling, Polyhedron &poly)
		{
			const int n_values = poly.n_nodes * poly.n_dims;
			for(int i = 0; i < n_values; ++i) {
				poly.points[i] *= scaling;
			}
		}

		template<int Dimensions>
		bool Assemble(
			const std::shared_ptr<MeshBase> &master,
			const std::shared_ptr<MeshBase> &slave,
			const std::shared_ptr<DofMap> &dof_master,
			const std::shared_ptr<DofMap> &dof_slave,
			const unsigned int &from_var_num,
			const unsigned int &to_var_num,
			DSMatrixd &B,
			const moonolith::SearchSettings &settings,
			bool  use_biorth,
			const std::vector< std::pair<int, int> > &tags,
			int n_var)
		{
			const bool use_interpolation = false;

			moonolith::Communicator comm(master->comm().get());

			const int var_num_slave = to_var_num;
			std::shared_ptr<FESpacesAdapter> local_fun_spaces = std::make_shared<FESpacesAdapter>(master, slave, dof_master, dof_slave,from_var_num,to_var_num);

			Local2Global local2global(use_interpolation);

			std::shared_ptr<LocalAssembler> assembler;
			if(use_interpolation) {
				// assembler = std::make_shared<InterpolationLocalAssembler>(master->mesh_dimension(), use_biorth);
			} else {
				assembler = std::make_shared<L2LocalAssembler>(master->mesh_dimension(), use_biorth);
			}


			libMesh::DenseMatrix<libMesh::Real> elemmat;
			libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;

			libMesh::Real local_element_matrices_sum = 0.0;

			moonolith::SparseMatrix<double> mat_buffer(comm);
			mat_buffer.set_size(dof_slave->n_dofs(), dof_master->n_dofs());

			auto fun = [&](
				const VElementAdapter<Dimensions> &master,
				const VElementAdapter<Dimensions> &slave) -> bool {

				const auto &master_mesh  = master.space();;

				const auto &slave_mesh = slave.space();

				const int src_index  = master.element();

				const int dest_index = slave.element();

				auto &master_el  = *master_mesh.elem(src_index);

				auto &slave_el = *slave_mesh.elem(dest_index);

				const int dim = master_mesh.mesh_dimension();

				bool vol2surf = false;

				auto master_type = dof_master->variable(0).type();
				auto slave_type  = dof_slave->variable(0).type();

				elemmat.zero();
				if(assembler->assemble(master_el, master_type, slave_el, slave_type, elemmat)) {
					auto partial_sum = std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(), libMesh::Real(0.0));
					local_element_matrices_sum += partial_sum;

					const auto &master_dofs = master.dof_map();
					const auto &slave_dofs  = slave.dof_map();

					local2global.apply(master_dofs, slave_dofs, elemmat, mat_buffer);
					return true;
				} else {
					return false;
				}
			};

			if(!Assemble<Dimensions>(master, slave, dof_master, dof_slave, from_var_num, to_var_num, fun, settings, use_biorth, tags, n_var)) {
				return false;
			}

			double total_intersection_volume = 0.;
			{
				auto l2_assembler = std::dynamic_pointer_cast<L2LocalAssembler>(assembler);
				if(l2_assembler) {
					total_intersection_volume = l2_assembler->get_q_builder().get_total_intersection_volume();
				}
			}

			double volumes[2] = { local_element_matrices_sum,  total_intersection_volume };

			comm.all_reduce(volumes, 2, moonolith::MPISum());

			const processor_id_type master_proc_id  = master->processor_id();

			const dof_id_type n_dofs_on_proc_master = dof_master->n_local_dofs();

			const processor_id_type slave_proc_id   = slave->processor_id();

			const dof_id_type n_dofs_on_proc_slave  =dof_slave->n_local_dofs();

			const int n_dofs_on_proc_print  = dof_slave->n_local_dofs();

			if(comm.is_root()) {
				std::cout << "sum(B): " << volumes[0] << ", vol(I): " << volumes[1] << std::endl;
			}

			std::vector<moonolith::Integer> ownershipRangesMaster(comm.size()+1, 0);
			std::vector<moonolith::Integer> ownershipRangesSlave(comm.size()+1, 0);

			ownershipRangesMaster[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_master);
			ownershipRangesSlave[comm.rank()+1]  += static_cast<unsigned int>(n_dofs_on_proc_slave);

			comm.all_reduce(&ownershipRangesMaster[0], ownershipRangesMaster.size(), moonolith::MPISum());
			comm.all_reduce(&ownershipRangesSlave[0],  ownershipRangesSlave.size(),  moonolith::MPISum());

			std::partial_sum(ownershipRangesMaster.begin(), ownershipRangesMaster.end(), ownershipRangesMaster.begin());
			std::partial_sum(ownershipRangesSlave.begin(), ownershipRangesSlave.end(), ownershipRangesSlave.begin());


			int dim = master->mesh_dimension();

			moonolith::Redistribute< moonolith::SparseMatrix<double> > redist(comm.get_mpi_comm());

			redist.apply(ownershipRangesSlave, mat_buffer, moonolith::AddAssign<double>());

			assert(ownershipRangesSlave.empty() == ownershipRangesMaster.empty() || ownershipRangesMaster.empty());

			SizeType  mMaxRowEntries = mat_buffer.local_max_entries_x_col();

			comm.all_reduce(&mMaxRowEntries, 1, moonolith::MPIMax());

			const SizeType local_range_slave_range  = ownershipRangesSlave [comm.rank()+1] - ownershipRangesSlave [comm.rank()];
			const SizeType local_range_master_range = ownershipRangesMaster[comm.rank()+1] - ownershipRangesMaster[comm.rank()];

			DSMatrixd B_x = utopia::local_sparse(local_range_slave_range, local_range_master_range, mMaxRowEntries);

			{
				utopia::Write<utopia::DSMatrixd> write(B_x);
				for (auto it = mat_buffer.iter(); it; ++it) {
					B_x.set(it.row(), it.col(), *it);

				}
			}

			auto s_B_x = local_size(B_x);

			B = local_sparse(s_B_x.get(0), s_B_x.get(1), n_var * mMaxRowEntries);

			utopia::Write<DSMatrixd> w_B(B);
			utopia::each_read(B_x, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
				for(utopia::SizeType d = 0; d < n_var; ++d) {
					B.set(i+d, j+d, value);
				}
			});

			return true;
		}

		bool assemble_volume_transfer(
			moonolith::Communicator &comm,
			const std::shared_ptr<MeshBase> &master,
			const std::shared_ptr<MeshBase> &slave,
			const std::shared_ptr<DofMap> &dof_master,
			const std::shared_ptr<DofMap> &dof_slave,
			const unsigned int &from_var_num,
			const unsigned int &to_var_num,
			bool use_biorth,
			int n_var, DSMatrixd &B,
			const std::vector< std::pair<int, int> > &tags)
		{
		///////////////////////////
			if(Utopia::instance().verbose()) {
				moonolith::root_describe("---------------------------------------\n"
					"begin: utopia::assemble_volume_transfer",
					comm, std::cout);
			}

			Chrono c;
			c.start();
		///////////////////////////

			moonolith::SearchSettings settings;
		// settings.verbosity_level = 3;

			bool ok = false;
			if(master->mesh_dimension() == 2) {
				ok = utopia::Assemble<2>(master, slave, dof_master, dof_slave, from_var_num,  to_var_num, B, settings,use_biorth, tags, n_var);
			} else if(master->mesh_dimension() == 3) {
				ok = utopia::Assemble<3>(master, slave, dof_master, dof_slave, from_var_num,  to_var_num, B, settings,use_biorth, tags, n_var);
			} else {
				assert(false && "Dimension not supported!");
			}

		///////////////////////////
			c.stop();

			if(Utopia::instance().verbose()) {
				std::stringstream ss;
				ss << "end: utopia::assemble_volume_transfer\n";
				ss << c;
				ss << "---------------------------------------";
				moonolith::root_describe(ss.str(), comm, std::cout);
			}
		///////////////////////////

			return ok;
		}
	}




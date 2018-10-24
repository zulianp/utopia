#ifndef UTOPIA_GRID_MESH_TRANSFER_HPP
#define UTOPIA_GRID_MESH_TRANSFER_HPP

#include "utopia.hpp"
#include "utopia_Grid.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_QMortarBuilder.hpp"
#include "MortarAssemble.hpp"

#include "utopia_intersector.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_Socket.hpp"

#include "libmesh/serial_mesh.h"
#include <functional>
#include <array>
#include <vector>

namespace utopia {

	
	class Grid2MeshTransferAssembler {
	public:
		using FunctionSpace = utopia::LibMeshFunctionSpace;
		using SparseMatrix  = utopia::USparseMatrix;
		using MeshBase      = libMesh::MeshBase;
		using DofMap        = libMesh::DofMap;

		using ElementMatrix = LocalAssembler::Matrix;

		Grid2MeshTransferAssembler(
			const std::shared_ptr<LocalAssembler> &assembler,
			const std::shared_ptr<Local2Global>   &local2global)
		: assembler_(assembler), local2global_(local2global)
		{}

		template<int Dim>
		bool assemble(
			const Grid<Dim> &from_mesh,
			const std::vector<long> &ownership_ranges,
			const std::shared_ptr<MeshBase> &to_mesh,
			const std::shared_ptr<DofMap>   &to_dofs,
			std::vector<std::shared_ptr<SparseMatrix>> &mats,
			const TransferOptions &opts = TransferOptions())
		{
			using Vector  = typename Grid<Dim>::Vector;
			using Array   = typename Grid<Dim>::Array;
			using Index   = typename Grid<Dim>::Index;
			using Scalar  = typename Grid<Dim>::Scalar;
			using Integer = typename Grid<Dim>::Integer;

			moonolith::Communicator comm(to_mesh->comm().get());

			mats.resize(assembler_->n_forms());

			pre_assemble(comm, from_mesh.n_nodes(), to_dofs->n_dofs());
			Index index;
			Vector emin, emax;
			std::vector<libMesh::dof_id_type> temp_slave_dofs;
			std::vector<long> master_dofs, slave_dofs;

			auto from_elem_type = grid_elem_type(from_mesh);

			for(auto it = to_mesh->active_local_elements_begin(); it != to_mesh->active_local_elements_end(); ++it) {
				const auto &e = **it;

				for(int i = 0; i < Dim; ++i) {
					emin[i] = std::numeric_limits<Scalar>::max();
					emax[i] = -std::numeric_limits<Scalar>::max();
				}

				for(auto k = 0; k < e.n_nodes(); ++k) {
					for(int i = 0; i < Dim; ++i) {
						emin[i] = std::min(emin[i], e.node_ref(k)(i));
						emax[i] = std::max(emax[i], e.node_ref(k)(i));
					}
				}

				from_mesh.elements_in_range(emin, emax, index);
				if(index.empty()) continue;

				for(auto ind : index) {
					// polyhedron(grid, ind, from_poly);
					auto temp_mesh = build_element_mesh(to_mesh->comm(), from_mesh, ind);
					auto grid_elem = temp_mesh->elem(0);

					for(auto &mat_i : elemmat) {
						mat_i.zero();
					}

					// temp_mesh->prepare_for_use();
					// plot_mesh(*temp_mesh, "grid/m" + std::to_string(ind));

					if(assembler_->assemble(
						*grid_elem,
						from_elem_type,
						e,
						to_dofs->variable_type(opts.to_var_num),
						elemmat)) {

						++n_intersections_;
					
						for(std::size_t i = 0; i < elemmat.size(); ++i) {	
							auto &mat_i = elemmat[i];
							auto partial_sum = std::accumulate(mat_i.get_values().begin(), mat_i.get_values().end(), libMesh::Real(0.0));

							assert(!std::isnan(partial_sum));
							local_element_matrices_sum[i] += partial_sum;

							switch(assembler_->type(i)) {
								
								case LocalAssembler::MASTER_X_SLAVE: 
								{
									from_mesh.dofs(ind, master_dofs);
									to_dofs->dof_indices(&e, temp_slave_dofs);
									slave_dofs.clear(); 
									slave_dofs.insert(slave_dofs.begin(), temp_slave_dofs.begin(), temp_slave_dofs.end());

									local2global_->apply(master_dofs, slave_dofs, elemmat[i], *mat_buffer[i]);
									break;
								}

								case LocalAssembler::SLAVE_X_SLAVE:
								{
									to_dofs->dof_indices(&e, temp_slave_dofs);
									slave_dofs.clear(); 
									slave_dofs.insert(slave_dofs.begin(), temp_slave_dofs.begin(), temp_slave_dofs.end());

									local2global_->apply(slave_dofs, slave_dofs, elemmat[i], *mat_buffer[i]);
									break;
								}

								case LocalAssembler::MASTER_X_MASTER:
								{
									from_mesh.dofs(ind, master_dofs);

									local2global_->apply(master_dofs, master_dofs, elemmat[i], *mat_buffer[i]);
									break;
								}

								case LocalAssembler::SLAVE_X_MASTER: 
								{
									from_mesh.dofs(ind, master_dofs);
									to_dofs->dof_indices(&e, temp_slave_dofs);
									slave_dofs.clear(); 
									slave_dofs.insert(slave_dofs.begin(), temp_slave_dofs.begin(), temp_slave_dofs.end());

									local2global_->apply(slave_dofs, master_dofs, elemmat[i], *mat_buffer[i]);
									break;
								}

								default:
								{
									assert(false);
									break;
								}
							}
						}


					}  else {
						
				
					}
				}
			}

			auto n_local_dofs_from = ownership_ranges[comm.rank() + 1] - ownership_ranges[comm.rank()];

			for(std::size_t i = 0; i < mat_buffer.size(); ++i) {
				post_assemble(comm, n_local_dofs_from, to_dofs, opts, mats, i);
			}

			print_stats(comm);
			return true;
		}

		void pre_assemble(moonolith::Communicator &comm, const std::size_t from_n_dofs, const std::size_t to_n_dofs)
		{
			init_buffers(assembler_->n_forms());
			const std::size_t n_forms = mat_buffer.size();

			for(std::size_t i = 0; i < n_forms; ++i) {
				mat_buffer[i] = std::make_shared< moonolith::SparseMatrix<double> >(comm);
				local_element_matrices_sum[i] = 0.;

				switch(assembler_->type(i)) {
					
					case LocalAssembler::MASTER_X_SLAVE: 
					{
						mat_buffer[i]->set_size(to_n_dofs, from_n_dofs);
						break;
					}

					case LocalAssembler::SLAVE_X_SLAVE:
					{
						mat_buffer[i]->set_size(to_n_dofs, to_n_dofs);
						break;
					}

					case LocalAssembler::MASTER_X_MASTER:
					{
						mat_buffer[i]->set_size(from_n_dofs, from_n_dofs);
						break;
					}

					case LocalAssembler::SLAVE_X_MASTER: 
					{
						mat_buffer[i]->set_size(from_n_dofs, to_n_dofs);
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

		void init_buffers(const SizeType n)
		{
			assert(n > 0);
			n_intersections_ = 0;
			mat_buffer.resize(n);
			elemmat.resize(n);
			local_element_matrices_sum.resize(n);
		}

		void post_assemble(
			moonolith::Communicator &comm,
			const std::size_t n_local_dofs_from,
			const std::shared_ptr<DofMap> &to_dofs,
			const TransferOptions &opts,
			std::vector<std::shared_ptr<SparseMatrix>> &mats,
			std::size_t buffer_num)
		{
			if(!mats[buffer_num]) {
				mats[buffer_num] = std::make_shared<SparseMatrix>();
			}

			SparseMatrix &mat = *mats[buffer_num];

			libMesh::dof_id_type n_dofs_on_proc_trial = 0;
			libMesh::dof_id_type n_dofs_on_proc_test  = 0;

			switch(assembler_->type(buffer_num)) {
				
				case LocalAssembler::MASTER_X_SLAVE: 
				{
					n_dofs_on_proc_trial = n_local_dofs_from;
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
					n_dofs_on_proc_trial = n_local_dofs_from;
					n_dofs_on_proc_test  = n_local_dofs_from;
					break;
				}

				case LocalAssembler::SLAVE_X_MASTER: 
				{
					n_dofs_on_proc_trial = to_dofs->n_local_dofs();
					n_dofs_on_proc_test  = n_local_dofs_from;
					break;
				}

				default:
				{
					assert(false);
					break;
				}
			}

			local2global_->redistribute(comm, n_dofs_on_proc_trial, n_dofs_on_proc_test, *mat_buffer[buffer_num]);

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

		template<int Dim>
		inline libMesh::FEType grid_elem_type(const Grid<Dim> &grid)
		{
			libMesh::FEType type;
			// type.family = libMesh::LAGRANGE;
			// type.order = libMesh::FIRST;
			return type;
		}

		inline std::unique_ptr<libMesh::SerialMesh> build_element_mesh(
			const libMesh::Parallel::Communicator &comm,
			const Grid<2> &grid,
			Grid<2>::Integer cell_index)
		{
			auto index_min = grid.element_index(cell_index);
			auto index_max = index_min;

			for(auto &i : index_max) {
				++i;
			}

			auto p_min = grid.point(index_min);
			auto p_max = grid.point(index_max);

			auto mesh = make_unique<libMesh::SerialMesh>(comm, 2);
			mesh->reserve_nodes(4);
			
			libMesh::Point p;
			//point 0
			p(0) = p_min.x;
			p(1) = p_min.y;
			mesh->add_point(p);

			//point 1
			p(0) = p_max.x;
			p(1) = p_min.y;
			mesh->add_point(p);

			//point 2
			p(0) = p_max.x;
			p(1) = p_max.y;
			mesh->add_point(p);

			//point 3
			p(0) = p_min.x;
			p(1) = p_max.y;
			mesh->add_point(p);

			auto elem = libMesh::Elem::build(libMesh::QUAD4);

			for (int i = 0; i < 4; ++i) {
				elem->set_node(i) = & mesh->node(i);
			}


			mesh->add_elem(elem.release());
			return std::move(mesh);
		}

		inline std::unique_ptr<libMesh::SerialMesh> build_element_mesh(
			const libMesh::Parallel::Communicator &comm,
			const Grid<3> &grid,
			Grid<3>::Integer cell_index)
		{
			
			auto index_min = grid.element_index(cell_index);
			auto index_max = index_min;

			for(auto &i : index_max) {
				++i;
			}

			auto p_min = grid.point(index_min);
			auto p_max = grid.point(index_max);

				/*
				* HEX8:   7        6
				*         o--------o
				*        /:       /|
				*       / :      / |
				*    4 /  :   5 /  |
				*     o--------o   |
				*     |   o....|...o 2
				*     |  .3    |  /
				*     | .      | /
				*     |.       |/
				*     o--------o
				*     0        1
				*/

			auto mesh = make_unique<libMesh::SerialMesh>(comm, 3);
			mesh->reserve_nodes(8);
			
			libMesh::Point p;
			//point 0
			p(0) = p_min.x;
			p(1) = p_min.y;
			p(2) = p_min.z;
			mesh->add_point(p);

			//point 1
			p(0) = p_max.x;
			p(1) = p_min.y;
			p(2) = p_min.z;
			mesh->add_point(p);

			//point 2
			p(0) = p_max.x;
			p(1) = p_max.y;
			p(2) = p_min.z;
			mesh->add_point(p);

			//point 3
			p(0) = p_min.x;
			p(1) = p_max.y;
			p(2) = p_min.z;
			mesh->add_point(p);

			//point 4
			p(0) = p_min.x;
			p(1) = p_min.y;
			p(2) = p_max.z;
			mesh->add_point(p);

			//point 5
			p(0) = p_max.x;
			p(1) = p_min.y;
			p(2) = p_max.z;
			mesh->add_point(p);

			//point 6
			p(0) = p_max.x;
			p(1) = p_max.y;
			p(2) = p_max.z;
			mesh->add_point(p);

			//point 7
			p(0) = p_min.x;
			p(1) = p_max.y;
			p(2) = p_max.z;
			mesh->add_point(p);

			auto elem = libMesh::Elem::build(libMesh::HEX8);

			for (int i = 0; i < 8; ++i) {
				elem->set_node(i) = & mesh->node(i);
			}


			mesh->add_elem(elem.release());
			return std::move(mesh);
		}
		
		void print_stats(moonolith::Communicator &comm)
		{
			{
				auto l2_assembler = std::dynamic_pointer_cast<L2LocalAssembler>(assembler_);
				if(l2_assembler) {
					double total_intersection_volume = l2_assembler->get_q_builder().get_total_intersection_volume();

					double volumes[2] = { local_element_matrices_sum[0], total_intersection_volume };
					comm.all_reduce(volumes, 2, moonolith::MPISum());
					comm.all_reduce(&n_intersections_, 1, moonolith::MPISum());

					if(comm.is_root()) {
						std::cout << "sum(B): " 
								  << volumes[0] 
								  << ", vol(I): " 
								  << volumes[1]
								  << ", n_intersections: "
								  << n_intersections_ 
								  << std::endl;
					}

					// moonolith::root_describe("vol(I_local) : " + std::to_string(total_intersection_volume), comm, std::cout);
				}
			}
		}

	public:
		std::shared_ptr<LocalAssembler> assembler_;
		std::shared_ptr<Local2Global> local2global_;
		
		std::vector< libMesh::DenseMatrix<libMesh::Real> > elemmat;
		std::vector< libMesh::Real > local_element_matrices_sum;

		std::vector< std::shared_ptr< moonolith::SparseMatrix<double> > > mat_buffer;
		long n_intersections_;
	};
}

#endif //UTOPIA_GRID_MESH_TRANSFER_HPP


#ifndef UTOPIA_GRID_MESH_TRANSFER_HPP
#define UTOPIA_GRID_MESH_TRANSFER_HPP


#include "utopia.hpp"
#include "utopia_TransferAssembler.hpp"
#include "MortarAssemble.hpp"

#include "utopia_intersector.hpp"
#include "moonolith_static_math.hpp"
#include <functional>
#include <array>
#include <vector>

namespace utopia {

	//grid defined in unit-cube
	template<int Dim>
	class Grid {
	public:
		using Scalar  = double;
		using Integer = long;
		using SizeT   = std::size_t;

		using Vector  = moonolith::Vector<Scalar, Dim>;
		using Array   = std::array<Integer, Dim>;
		using Mapping = std::function<Vector (const Vector &)>;
		using Index   = std::vector<Integer>;

		inline Array index(const Integer hash) const
		{
			assert(false && "implement me");
			return Array();
		}

		inline Vector point(const Array &index) const
		{
			Vector ret;
			for(int i = 0; i < Dim; ++i) {
				ret[i] = static_cast<Scalar>(index[i])/dims[i];
			}

			return map(ret);
		}

		inline Integer hash(const Vector &y) const
		{
			Vector x = inverse_map(y);
			Integer result = floor(x[0] * dims[0]);

			Integer total_dim = dims[0];

			for(int i = 1; i < Dim; ++i) {
				result *= dims[i];
				result += floor(x[i] * dims[i]);
				total_dim *= dims[i];
			}

			if(result >= total_dim || result < 0) {
				printf("error -> %ld\n", result);
			}

			return result;
		}

		inline Integer hash(const Array &coord) const
		{
			return hash_from_index(coord);
		}

		inline Integer hash_from_index(const Array &coord) const
		{
			Integer result   = coord[0];
			Integer total_dim = dims[0];

			for(int i = 1; i < Dim; ++i) {
				result *= dims[i];
				result += coord[i];
				total_dim *= dims[i];
			}

			return result;
		}

		template<typename IntT>
		void dofs(const Integer cell_index, std::vector<IntT> &dof_indices) const
		{
			dof_indices.resize(moonolith::StaticPow<2, Dim>::value);

			auto imin = index(cell_index);
			auto imax = imin;

			for(auto &i : imax) {++i; }

			if(Dim == 3) {

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

				dof_indices[0] = hash_from_index(imin);
				dof_indices[1] = hash_from_index({ imax[0], imin[1], imin[2] });
				dof_indices[2] = hash_from_index({ imax[0], imax[1], imin[2] });
				dof_indices[3] = hash_from_index({ imin[0], imax[1], imin[2] });
				dof_indices[4] = hash_from_index({ imin[0], imin[1], imax[2] });
				dof_indices[5] = hash_from_index({ imax[0], imin[1], imax[2] });
				dof_indices[6] = hash_from_index(imax);
				dof_indices[7] = hash_from_index({ imin[0], imax[1], imax[2] });
			} else {
				assert(false && "implement me");
			}
		}

		

		void hash_range(
			const Vector &min,
			const Vector &max, Index &hashes) const
		{
			Array imin, imax, offsets;

			auto im_min = inverse_map(min);
			auto im_max = inverse_map(max);

			//generate tensor indices
			for(int i = 0; i < Dim; ++i) {
				imin[i] = floor(im_min[i] * dims[i]);
			}

			for(int i = 0; i < Dim; ++i) {
				imax[i] = floor(im_max[i] * dims[i]);
			}

			for(int i = 0; i < Dim; ++i) {
				offsets[i] = imax[i] - imin[i];
			}

			//FIXME make more general for greater dim and extract to class
			Array coord;
			if(Dim == 1) {
				for(int i = imin[0]; i <= imax[0]; ++i) {
					coord[0] = i;
					hashes.push_back(hash(coord));
				}

			} else if(Dim == 2) {
				for(int i = imin[0]; i <= imax[0]; ++i) {
					for(int j = imin[1]; j <= imax[1]; ++j) {
						coord[0] = i;
						coord[1] = j;
						hashes.push_back(hash(coord));
					}
				}
			} else if(Dim == 3) {
				for(int i = imin[0]; i <= imax[0]; ++i) {
					for(int j = imin[1]; j <= imax[1]; ++j) {
						for(int k = imin[2]; k <= imax[2]; ++k) {
							coord[0] = i;
							coord[1] = j;
							coord[2] = k;
							hashes.push_back(hash(coord));
						}
					}
				}
			} else {
				assert(false && "dim > 3 not supported yet!");
			}

			assert(!hashes.empty());
		}


		Grid()
		{
			map = [](const Vector &x) -> Vector { return x; };
			inverse_map = map;
			std::fill(std::begin(dims), std::end(dims), 0);
		}

		inline Integer n_nodes() const
		{
			Integer ret = 1;

			for(auto d : dims) {
				ret *= d;
			}

			return ret;
		}

		Array dims;

		Mapping map;
		Mapping inverse_map;
	};

	template<int Dim>
	class GridMeshTransfer {
	public:
		using Vector = typename Grid<Dim>::Vector;
		using Array  = typename Grid<Dim>::Array;
		using Index  = typename Grid<Dim>::Index;
		using Scalar = typename Grid<Dim>::Scalar;
		using Integer = typename Grid<Dim>::Integer;

		using FunctionSpace = utopia::LibMeshFunctionSpace;
		using SparseMatrix  = utopia::USparseMatrix;
		using MeshBase      = libMesh::MeshBase;
		using DofMap        = libMesh::DofMap;

		using ElementMatrix = LocalAssembler::Matrix;

		GridMeshTransfer(
			const std::shared_ptr<LocalAssembler> &assembler,
			const std::shared_ptr<Local2Global>   &local2global)
		: assembler_(assembler), local2global_(local2global)
		{

		}

		bool assemble(
			const Grid<Dim> &from_mesh,
			const std::vector<long> &ownership_ranges,
			const std::shared_ptr<MeshBase> &to_mesh,
			const std::shared_ptr<DofMap>   &to_dofs,
			std::vector<std::shared_ptr<SparseMatrix> > &B,
			const TransferOptions &opts = TransferOptions())
		{
			moonolith::Communicator comm(to_mesh->comm().get());

			pre_assemble(comm, from_mesh.n_nodes(), to_dofs->n_dofs());

			std::vector<libMesh::Node> nodes;

			Index index;
			Vector emin, emax;
			Polyhedron from_poly, to_poly;

			// std::vector<ElementMatrix> elemmat;
			// std::vector<double> local_element_matrices_sum(assembler_->n_forms(), 0.);

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

				from_mesh.hash_range(emin, emax, index);
				if(index.empty()) continue;

				for(auto ind : index) {
					// polyhedron(grid, ind, from_poly);
					auto grid_elem = build_element(from_mesh, ind, nodes);

					for(auto &mat_i : elemmat) {
						mat_i.zero();
					}

					if(assembler_->assemble(
						*grid_elem,
						from_elem_type,
						e,
						to_dofs->variable_type(opts.to_var_num),
						elemmat)) {

						for(std::size_t i = 0; i < elemmat.size(); ++i) {	
							auto &mat_i = elemmat[i];
							auto partial_sum = std::accumulate(mat_i.get_values().begin(), mat_i.get_values().end(), libMesh::Real(0.0));
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
				post_assemble(comm, n_local_dofs_from, to_dofs, opts, i);
			}

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
			mat_buffer.resize(n);
			elemmat.resize(n);
			local_element_matrices_sum.resize(n);
		}

		void post_assemble(
			moonolith::Communicator &comm,
			const std::size_t n_local_dofs_from,
			const std::shared_ptr<DofMap> &to_dofs,
			const TransferOptions &opts,
			std::size_t buffer_num)
		{
			SparseMatrix &mat = *mats_[buffer_num];

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


		inline libMesh::FEType grid_elem_type(const Grid<3> &grid)
		{
			libMesh::FEType type;
			type.family = libMesh::LAGRANGE;
			type.order = libMesh::FIRST;
			return type;
		}

		inline libMesh::UniquePtr<libMesh::Elem> build_element(
			const Grid<3> &grid,
			Integer cell_index,
			std::vector<libMesh::Node> &nodes)
		{
			auto elem = libMesh::Elem::build(libMesh::HEX8);

			Array index_min = grid.index(cell_index);
			Array index_max = index_min;

			for(auto &i : index_max) {
				++i;
			}

			Vector p_min = grid.point(index_min);
			Vector p_max = grid.point(index_max);

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

			nodes.resize(8);
			
			//point 0
			nodes[0](0) = p_min.x;
			nodes[0](0) = p_min.y;
			nodes[0](0) = p_min.z;

			//point 1
			nodes[1](0) = p_max.x;
			nodes[1](1) = p_min.y;
			nodes[1](2) = p_min.z;

			//point 2
			nodes[2](0) = p_max.x;
			nodes[2](1) = p_max.y;
			nodes[2](2) = p_min.z;

			//point 3
			nodes[3](0) = p_min.x;
			nodes[3](1) = p_min.y;
			nodes[3](2) = p_max.z;

			//point 4
			nodes[4](0) = p_min.x;
			nodes[4](1) = p_max.y;
			nodes[4](2) = p_min.z;

			//point 5
			nodes[5](0) = p_max.x;
			nodes[5](1) = p_min.y;
			nodes[5](2) = p_max.z;

			//point 6
			nodes[6](0) = p_max.x;
			nodes[6](1) = p_max.y;
			nodes[6](2) = p_max.z;

			//point 7
			nodes[7](0) = p_min.x;
			nodes[7](1) = p_max.y;
			nodes[7](2) = p_max.z;

			for (int i = 0; i < 8; ++i) {
				elem->set_node(i) = &nodes[i];
			}

			return std::move(elem);

		}

			//compatibility method
		inline void polyhedron(
			const Grid<Dim> &grid,
			const Integer cell_index,
			Polyhedron &poly)
		{
			assert(Dim == 3);

			poly.n_elements = 6;
			poly.n_dims  = Dim;
			poly.n_nodes = 8;
			poly.type    = P_MESH_TYPE_HEX;

			Array index_min = grid.index(cell_index);
			Array index_max = index_min;

			for(auto &i : index_max) {
				++i;
			}

			Vector p_min = grid.point(index_min);
			Vector p_max = grid.point(index_max);

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

				//point 0
			poly.points[0] = p_min.x;
			poly.points[1] = p_min.y;
			poly.points[2] = p_min.z;

				//point 1
			poly.points[3] = p_max.x;
			poly.points[4] = p_min.y;
			poly.points[5] = p_min.z;

				//point 2
			poly.points[6] = p_max.x;
			poly.points[7] = p_max.y;
			poly.points[8] = p_min.z;

				//point 3
			poly.points[6] = p_min.x;
			poly.points[7] = p_min.y;
			poly.points[8] = p_max.z;

				//point 4
			poly.points[6] = p_min.x;
			poly.points[7] = p_max.y;
			poly.points[8] = p_min.z;

				//point 5
			poly.points[9]  = p_max.x;
			poly.points[10] = p_min.y;
			poly.points[11] = p_max.z;

				//point 6
			poly.points[12] = p_max.x;
			poly.points[13] = p_max.y;
			poly.points[14] = p_max.z;

				//point 7
			poly.points[15] = p_min.x;
			poly.points[16] = p_max.y;
			poly.points[17] = p_max.z;

				//el pointer
			poly.el_ptr[0] = 0;
			poly.el_ptr[1] = 4;
			poly.el_ptr[2] = 8;
			poly.el_ptr[3] = 12;
			poly.el_ptr[4] = 16;
			poly.el_ptr[5] = 20;
			poly.el_ptr[6] = 24;

				//face 0
			poly.el_index[0] = 0;
			poly.el_index[1] = 1;
			poly.el_index[2] = 5;
			poly.el_index[3] = 4;

				//face 1
			poly.el_index[4] = 1;
			poly.el_index[5] = 2;
			poly.el_index[6] = 6;
			poly.el_index[7] = 5;

				//face 2
			poly.el_index[8]  = 3;
			poly.el_index[9]  = 7;
			poly.el_index[10] = 6;
			poly.el_index[11] = 2;

				//face 3
			poly.el_index[12] = 0;
			poly.el_index[13] = 4;
			poly.el_index[14] = 7;
			poly.el_index[15] = 3;

				//face 4
			poly.el_index[16] = 2;
			poly.el_index[17] = 1;
			poly.el_index[18] = 0;
			poly.el_index[19] = 3;

				//face 5
			poly.el_index[20] = 6;
			poly.el_index[21] = 7;
			poly.el_index[22] = 4;
			poly.el_index[23] = 5;

			poly.type = P_MESH_TYPE_HEX;
		}

	public:
		std::shared_ptr<LocalAssembler> assembler_;
		std::shared_ptr<Local2Global> local2global_;
		
		std::vector< libMesh::DenseMatrix<libMesh::Real> > elemmat;
		std::vector< libMesh::Real > local_element_matrices_sum;

		std::vector< std::shared_ptr< moonolith::SparseMatrix<double> > > mat_buffer;
		std::vector<std::shared_ptr<SparseMatrix>> mats_;
	};
}

#endif //UTOPIA_GRID_MESH_TRANSFER_HPP


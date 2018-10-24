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

#include "utopia_private_PetrovGalerkinAssembler.hpp"

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
			auto n_local_dofs_from = ownership_ranges[comm.rank() + 1] - ownership_ranges[comm.rank()];

			pg_assembler_.initialize(
				comm,
				assembler_,
				local2global_,
				opts,
				from_mesh.n_nodes(),
				n_local_dofs_from,
				to_dofs->n_dofs(),
				to_dofs->n_local_dofs()
			);

			Index index;
			Vector emin, emax;
			std::vector<libMesh::dof_id_type> temp_slave_dofs;

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
					auto temp_mesh = build_element_mesh(to_mesh->comm(), from_mesh, ind);
					auto grid_elem = temp_mesh->elem(0);

					// temp_mesh->prepare_for_use();
					// plot_mesh(*temp_mesh, "grid/m" + std::to_string(ind));
					pg_assembler_.assemble(
						*grid_elem,
						from_elem_type,
						e,
						to_dofs->variable_type(opts.to_var_num),
						[&](std::vector<long> &master_dofs, std::vector<long> &slave_dofs) {
							from_mesh.dofs(ind, master_dofs);
							to_dofs->dof_indices(&e, temp_slave_dofs);
							slave_dofs.clear(); 
							slave_dofs.insert(slave_dofs.begin(), temp_slave_dofs.begin(), temp_slave_dofs.end());
						}
					);
				}
			}

			pg_assembler_.finalize(mats);
			pg_assembler_.print_stats();
			return true;
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

	public:
		std::shared_ptr<LocalAssembler> assembler_;
		std::shared_ptr<Local2Global> local2global_;
		private_::PetrovGalerkinAssembler pg_assembler_;


	};
}

#endif //UTOPIA_GRID_MESH_TRANSFER_HPP


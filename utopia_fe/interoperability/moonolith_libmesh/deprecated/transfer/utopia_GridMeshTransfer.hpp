#ifndef UTOPIA_GRID_MESH_TRANSFER_HPP
#define UTOPIA_GRID_MESH_TRANSFER_HPP

#include "MortarAssemble.hpp"
#include "utopia.hpp"
#include "utopia_Grid.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_QMortarBuilder.hpp"
#include "utopia_TransferAssembler.hpp"

#include "utopia_Socket.hpp"
#include "utopia_intersector.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_Voxel2Element.hpp"
#include "utopia_private_PetrovGalerkinAssembler.hpp"

#include <array>
#include <functional>
#include <vector>
#include "libmesh/serial_mesh.h"

namespace utopia {

    class Grid2MeshTransferAssembler {
    public:
        using FunctionSpace = utopia::LibMeshFunctionSpace;
        using SparseMatrix = utopia::USparseMatrix;
        using MeshBase = libMesh::MeshBase;
        using DofMap = libMesh::DofMap;

        using ElementMatrix = LocalAssembler::Matrix;

        Grid2MeshTransferAssembler(const std::shared_ptr<LocalAssembler> &assembler,
                                   const std::shared_ptr<Local2Global> &local2global)
            : assembler_(assembler), local2global_(local2global) {}

        template <int Dim>
        bool assemble(const Grid<Dim> &from_mesh,
                      const std::vector<long> &ownership_ranges,
                      const std::shared_ptr<MeshBase> &to_mesh,
                      const std::shared_ptr<DofMap> &to_dofs,
                      std::vector<std::shared_ptr<SparseMatrix>> &mats,
                      const TransferOptions &opts = TransferOptions()) {
            using Vector = typename Grid<Dim>::Vector;
            using Array = typename Grid<Dim>::Array;
            using Index = typename Grid<Dim>::Index;
            using Scalar = typename Grid<Dim>::Scalar;
            using Integer = typename Grid<Dim>::Integer;

            moonolith::Communicator comm(to_mesh->comm().get());

            if (Utopia::instance().verbose()) {
                moonolith::root_describe(
                    "---------------------------------------\n"
                    "begin: utopia::Grid2MeshTransferAssembler::assemble",
                    comm,
                    std::cout);
            }

            Chrono c;
            c.start();

            auto n_local_dofs_from = ownership_ranges[comm.rank() + 1] - ownership_ranges[comm.rank()];

            pg_assembler_.initialize(comm,
                                     assembler_,
                                     local2global_,
                                     opts,
                                     from_mesh.n_nodes(),
                                     n_local_dofs_from,
                                     to_dofs->n_dofs(),
                                     to_dofs->n_local_dofs());

            Index index;
            Vector emin, emax;
            std::vector<libMesh::dof_id_type> temp_slave_dofs;

            auto from_elem_type = grid_elem_type(from_mesh);

            for (auto it = to_mesh->active_local_elements_begin(); it != to_mesh->active_local_elements_end(); ++it) {
                const auto &e = **it;

                for (int i = 0; i < Dim; ++i) {
                    emin[i] = std::numeric_limits<Scalar>::max();
                    emax[i] = -std::numeric_limits<Scalar>::max();
                }

                for (auto k = 0; k < e.n_nodes(); ++k) {
                    for (int i = 0; i < Dim; ++i) {
                        emin[i] = std::min(emin[i], e.node_ref(k)(i));
                        emax[i] = std::max(emax[i], e.node_ref(k)(i));
                    }
                }

                from_mesh.elements_in_range(emin, emax, index);
                if (index.empty()) continue;

                for (auto ind : index) {
                    auto temp_mesh = Voxel2Element::build(to_mesh->comm(), from_mesh, ind);
                    auto grid_elem = utopia::elem_ptr(*temp_mesh, 0);

                    // temp_mesh->prepare_for_use();
                    // plot_mesh(*temp_mesh, "grid/m" + std::to_string(ind));
                    pg_assembler_.assemble(*grid_elem,
                                           from_elem_type,
                                           e,
                                           to_dofs->variable_type(opts.to_var_num),
                                           [&](std::vector<long> &master_dofs, std::vector<long> &slave_dofs) {
                                               from_mesh.dofs(ind, master_dofs);
                                               to_dofs->dof_indices(&e, temp_slave_dofs);
                                               slave_dofs.clear();
                                               slave_dofs.insert(
                                                   slave_dofs.begin(), temp_slave_dofs.begin(), temp_slave_dofs.end());
                                           });
                }
            }

            pg_assembler_.finalize(mats);
            pg_assembler_.print_stats();

            c.stop();

            if (Utopia::instance().verbose()) {
                std::stringstream ss;
                ss << "end: utopia::Grid2MeshTransferAssembler::assemble\n";
                ss << c;
                ss << "---------------------------------------";
                ss << "\n";
                assembler_->print_stats(ss);
                moonolith::root_describe(ss.str(), comm, std::cout);
            }
            return true;
        }

        template <int Dim>
        inline libMesh::FEType grid_elem_type(const Grid<Dim> &grid) {
            libMesh::FEType type;
            // type.family = libMesh::LAGRANGE;
            // type.order = libMesh::FIRST;
            return type;
        }

    public:
        std::shared_ptr<LocalAssembler> assembler_;
        std::shared_ptr<Local2Global> local2global_;
        private_::PetrovGalerkinAssembler pg_assembler_;
    };
}  // namespace utopia

#endif  // UTOPIA_GRID_MESH_TRANSFER_HPP

#ifndef UTOPIA_GRID_2_MESH_SURFACE_TRANSFER_ASSEMBLER_HPP
#define UTOPIA_GRID_2_MESH_SURFACE_TRANSFER_ASSEMBLER_HPP

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
#include <map>
#include <vector>
#include "libmesh/serial_mesh.h"

namespace utopia {

    class Grid2MeshSurfaceTransferAssembler {
    public:
        using FunctionSpace = utopia::LibMeshFunctionSpace;
        using SparseMatrix = utopia::USparseMatrix;
        using MeshBase = libMesh::MeshBase;
        using DofMap = libMesh::DofMap;

        using ElementMatrix = LocalAssembler::Matrix;

        template <int Dim>
        class Box {
        public:
            Box() {
                for (int d = 0; d < Dim; ++d) {
                    min[d] = std::numeric_limits<double>::max();
                    max[d] = -std::numeric_limits<double>::max();
                }
            }

            inline bool empty() const { return min[0] > max[0]; }

            typename Grid<Dim>::Vector min, max;
        };

        Grid2MeshSurfaceTransferAssembler(const std::shared_ptr<LocalAssembler> &assembler,
                                          const std::shared_ptr<Local2Global> &local2global)
            : assembler_(assembler), local2global_(local2global) {}

        static void side_dof_indices(const DofMap &dof_map,
                                     const int var_num,
                                     const libMesh::Elem &e,
                                     const int side_num,
                                     std::vector<long> &side_indices) {
            auto side_ptr = e.build_side_ptr(side_num);
            std::shared_ptr<Transform> trafo, side_trafo;

            if (e.dim() == 3) {
                trafo = std::make_shared<Transform3>(e);
                side_trafo = std::make_shared<Transform2>(*side_ptr);
            } else {
                assert(e.dim() == 2);
                trafo = std::make_shared<Transform2>(e);
                side_trafo = std::make_shared<Transform1>(*side_ptr);
            }

            QMortar q(e.dim());
            QMortar side_q(side_ptr->dim());

            q.get_points().resize(side_ptr->n_nodes());
            side_q.get_points().resize(side_ptr->n_nodes());

            for (auto i = 0; i < side_ptr->n_nodes(); ++i) {
                trafo->transform_to_reference(side_ptr->node_ref(i), q.get_points()[i]);
                side_trafo->transform_to_reference(side_ptr->node_ref(i), side_q.get_points()[i]);
            }

            auto fe = libMesh::FEBase::build(e.dim(), dof_map.variable_type(var_num));
            auto side_fe = libMesh::FEBase::build(side_ptr->dim(), dof_map.variable_type(var_num));

            auto &f = fe->get_phi();
            fe->attach_quadrature_rule(&q);
            fe->reinit(&e);

            auto &side_f = side_fe->get_phi();
            side_fe->attach_quadrature_rule(&side_q);
            side_fe->reinit(side_ptr.get());

            std::vector<libMesh::dof_id_type> elem_indices;
            dof_map.dof_indices(&e, elem_indices);

            side_indices.resize(side_f.size(), 0);

            auto n_qps = f[0].size();

            for (std::size_t i = 0; i < f.size(); ++i) {
                for (std::size_t j = 0; j < side_f.size(); ++j) {
                    double val = 0.;

                    for (std::size_t qp = 0; qp < n_qps; ++qp) {
                        val += f[i][qp] * side_f[j][qp];
                    }

                    if (std::abs(val) > 0.9) {
                        side_indices[j] = elem_indices[i];
                        break;
                    }
                }
            }
        }

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
                    "begin: utopia::Grid2MeshSurfaceTransferAssembler::assemble",
                    comm,
                    std::cout);
            }

            Chrono c;
            c.start();

            skipped_duplicate_intersections_ = 0;

            auto n_local_dofs_from = ownership_ranges[comm.rank() + 1] - ownership_ranges[comm.rank()];

            pg_assembler_.initialize(comm,
                                     assembler_,
                                     local2global_,
                                     opts,
                                     from_mesh.n_nodes(),
                                     n_local_dofs_from,
                                     to_dofs->n_dofs(),
                                     to_dofs->n_local_dofs());

            Array element_index;
            Index index;
            Vector emin, emax, erange;
            std::array<bool, Dim> is_planar;
            std::vector<std::set<Integer>> dual_graph;
            std::vector<Box<Dim>> boxes;
            std::vector<bool> is_coplanar;
            std::map<Integer, Integer> element_global_to_local_index;
            std::vector<bool> intersected;

            std::vector<libMesh::dof_id_type> temp_slave_dofs;

            auto from_elem_type = grid_elem_type(from_mesh);

            for (auto it = to_mesh->active_local_elements_begin(); it != to_mesh->active_local_elements_end(); ++it) {
                const auto &e = **it;

                for (uint side_num = 0; side_num < e.n_sides(); ++side_num) {
                    if (e.neighbor_ptr(side_num) != nullptr) {
                        continue;
                    }

                    auto side_elem_ptr = e.build_side_ptr(side_num);
                    auto &side_elem = *side_elem_ptr;

                    for (int i = 0; i < Dim; ++i) {
                        emin[i] = std::numeric_limits<Scalar>::max();
                        emax[i] = -std::numeric_limits<Scalar>::max();
                    }

                    for (auto k = 0; k < side_elem.n_nodes(); ++k) {
                        for (int i = 0; i < Dim; ++i) {
                            emin[i] = std::min(emin[i], side_elem.node_ref(k)(i));
                            emax[i] = std::max(emax[i], side_elem.node_ref(k)(i));
                        }
                    }

                    from_mesh.elements_in_range(emin, emax, index);
                    if (index.empty()) continue;

                    // detect face that is aligned or almost with side of voxel grid and mark voxel neighs
                    erange = emax - emin;

                    std::fill(std::begin(is_planar), std::end(is_planar), false);

                    int n_planar = 0;
                    double planar_tol = 1e-8;

                    for (int d = 0; d < Dim; ++d) {
                        if (approxeq(erange[d], 0., planar_tol)) {
                            ++n_planar;
                            is_planar[d] = true;
                        }
                    }

                    assert(n_planar <= 1);

                    if (n_planar > 0) {
                        // build local partial dual graph (only for planar dim)
                        // find adjacent elements
                        // check if there is a surface intersection and mark them
                        // check if surface intersection and assemble operator or discard duplicate

                        std::size_t n_candidates = index.size();
                        dual_graph.clear();
                        dual_graph.resize(n_candidates);

                        boxes.clear();
                        boxes.resize(n_candidates);

                        is_coplanar.clear();
                        is_coplanar.resize(n_candidates, false);

                        element_global_to_local_index.clear();

                        for (int d = 0; d < Dim; ++d) {
                            if (!is_planar[d]) continue;

                            for (std::size_t i = 0; i < n_candidates; ++i) {
                                auto ind = index[i];

                                if (boxes[i].empty()) {
                                    from_mesh.element_aabb(ind, boxes[i].min, boxes[i].max);
                                }

                                const bool is_min_coplanar = approxeq(emin[d], boxes[i].min[d], planar_tol);
                                const bool is_max_coplanar = approxeq(emin[d], boxes[i].max[d], planar_tol);

                                if (is_min_coplanar || is_max_coplanar) {
                                    // build dual graph for candidate i/ind for dimension d
                                    from_mesh.element_index(ind);
                                    element_index[d] += (is_min_coplanar) ? -1 : 1;

                                    if (from_mesh.element_is_valid(element_index)) {
                                        is_coplanar[i] = true;
                                        auto neigh = from_mesh.element_hash_from_index(element_index);
                                        element_global_to_local_index[ind] = i;
                                        dual_graph[i].insert(neigh);
                                    }
                                }
                            }
                        }

                        intersected.resize(n_candidates);
                        std::fill(std::begin(intersected), std::end(intersected), false);

                        for (std::size_t i = 0; i < n_candidates; ++i) {
                            auto ind = index[i];

                            bool skip_element = false;
                            if (is_coplanar[i]) {
                                for (auto neigh : dual_graph[i]) {
                                    auto it = element_global_to_local_index.find(neigh);
                                    assert(it != element_global_to_local_index.end());
                                    if (it == element_global_to_local_index.end()) continue;

                                    if (intersected[it->second]) {
                                        ++skipped_duplicate_intersections_;
                                        skip_element = true;
                                    }
                                }
                            }

                            if (skip_element) continue;

                            auto temp_mesh = Voxel2Element::build(to_mesh->comm(), from_mesh, ind);
                            auto grid_elem = temp_mesh->elem(0);

                            auto dof_fun = [&](std::vector<long> &master_dofs, std::vector<long> &slave_dofs) {
                                from_mesh.dofs(ind, master_dofs);

                                side_dof_indices(*to_dofs, opts.to_var_num, e, side_num, slave_dofs);
                            };

                            pg_assembler_.assemble(*grid_elem,
                                                   from_elem_type,
                                                   side_elem,
                                                   to_dofs->variable_type(opts.to_var_num),
                                                   dof_fun);
                        }

                    } else {
                        for (auto ind : index) {
                            auto temp_mesh = Voxel2Element::build(to_mesh->comm(), from_mesh, ind);
                            auto grid_elem = temp_mesh->elem(0);

                            // temp_mesh->prepare_for_use();
                            // plot_mesh(*temp_mesh, "grid/m" + std::to_string(ind));

                            auto dof_fun = [&](std::vector<long> &master_dofs, std::vector<long> &slave_dofs) {
                                from_mesh.dofs(ind, master_dofs);

                                side_dof_indices(*to_dofs, opts.to_var_num, e, side_num, slave_dofs);
                            };

                            pg_assembler_.assemble(*grid_elem,
                                                   from_elem_type,
                                                   side_elem,
                                                   to_dofs->variable_type(opts.to_var_num),
                                                   dof_fun);
                        }
                    }
                }
            }

            pg_assembler_.finalize(mats);
            pg_assembler_.print_stats();

            c.stop();

            if (Utopia::instance().verbose()) {
                std::stringstream ss;
                ss << "end: utopia::Grid2MeshSurfaceTransferAssembler::assemble\n";
                ss << c;
                ss << "\n";
                ss << "skipped duplicate intersections " << skipped_duplicate_intersections_ << "\n";

                assembler_->print_stats(ss);
                ss << "---------------------------------------";
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
        long skipped_duplicate_intersections_;
    };
}  // namespace utopia

#endif  // UTOPIA_GRID_2_MESH_SURFACE_TRANSFER_ASSEMBLER_HPP

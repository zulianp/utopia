#ifndef UTOPIA_KOKKOS_CUBE_HPP
#define UTOPIA_KOKKOS_CUBE_HPP

#include "utopia_Input.hpp"

#include <Kokkos_DynRankView.hpp>

namespace utopia {
    namespace kokkos {

        template <class Mesh>
        class Cube : public Configurable {
        public:
            using ExecutionSpace = typename Mesh::ExecutionSpace;
            using Range1 = Kokkos::RangePolicy<ExecutionSpace>;
            using Range3 = Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>;
            using Range2 = Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>;

            void read(Input &in) override {
                in.get("nx", n[0]);
                in.get("ny", n[1]);
                in.get("nz", n[2]);

                in.get("min_x", min[0]);
                in.get("min_y", min[1]);
                in.get("min_z", min[2]);

                in.get("max_x", max[0]);
                in.get("max_y", max[1]);
                in.get("max_z", max[2]);

                in.get("debug", debug);
            }

            template <class CrsGraph>
            void build_node_crs_graph(CrsGraph &crs_graph) {
                using RowPtr = typename CrsGraph::row_map_type::non_const_type;
                using ColIdx = typename CrsGraph::entries_type::non_const_type;
                using ColEntryType = typename CrsGraph::entries_type::value_type;

                int nx = n[0];
                int ny = n[1];
                int nz = n[2];

                int nrows = (nx + 1) * (ny + 1) * (nz + 1);

                RowPtr row_ptr("row_ptr", nrows + 1);
                ColIdx col_idx("col_idx", nrows * 27);

                Kokkos::parallel_for(
                    "Cube::build_node_crs_graph",
                    Range3({0, 0, 0}, {nx + 1, ny + 1, nz + 1}),
                    KOKKOS_LAMBDA(int i, int j, int k) {
                        int inode = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
                        int disp[3] = {-1, 0, 1};

                        int counter = 0;

                        // FIXME
                        for (int kk = 0; kk < 3; ++kk) {
                            int kkk = k + kk;
                            if (kkk < 0 || kkk > nz) continue;  // Avoid going outside

                            for (int jj = 0; jj < 3; ++jj) {
                                int jjj = j + jj;
                                if (jjj < 0 || jjj > ny) continue;  // Avoid going outside

                                for (int ii = 0; ii < 3; ++ii) {
                                    int iii = i + ii;
                                    if (iii < 0 || iii > nx) continue;  // Avoid going outside

                                    // int ineigh = (k + kk) * (ny + 1) * (nx + 1) + (j + jj) * (nx + 1) + (i + ii);
                                    counter++;
                                }
                            }
                        }

                        row_ptr(inode + 1) = counter;
                    });

                Kokkos::parallel_scan(
                    Range1(0, nrows + 1), KOKKOS_LAMBDA(const int &i, ColEntryType &upd, const bool &final) {
                        // Load old value in case we update it before accumulating
                        const ColEntryType val_i = row_ptr(i);

                        upd += val_i;

                        if (final) {
                            row_ptr(i) = upd;  // only update array on final pass
                        }
                    });

                Kokkos::parallel_for(
                    "Cube::build_node_crs_graph",
                    Range3({0, 0, 0}, {nx + 1, ny + 1, nz + 1}),
                    KOKKOS_LAMBDA(int i, int j, int k) {
                        int inode = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
                        int offset = row_ptr(inode);

                        int disp[3] = {-1, 0, 1};

                        int counter = 0;
                        for (int kk = 0; kk < 3; ++kk) {
                            int kkk = k + kk;
                            if (kkk < 0 || kkk > nz) continue;  // Avoid going outside

                            for (int jj = 0; jj < 3; ++jj) {
                                int jjj = j + jj;
                                if (jjj < 0 || jjj > ny) continue;  // Avoid going outside

                                for (int ii = 0; ii < 3; ++ii) {
                                    int iii = i + ii;
                                    if (iii < 0 || iii > nx) continue;  // Avoid going outside

                                    int ineigh = kkk * (ny + 1) * (nx + 1) + jjj * (nx + 1) + iii;
                                    col_idx(offset + counter) = ineigh;
                                    counter++;
                                }
                            }
                        }
                    });

                crs_graph = CrsGraph(col_idx, row_ptr);
            }

            void build(Mesh &mesh) {
                Chrono c;
                c.start();

                int nx = n[0];
                int ny = n[1];
                int nz = n[2];

                double min_x = min[0];
                double min_y = min[1];
                double min_z = min[2];

                double hx = max[0] - min_x;
                double hy = max[1] - min_y;
                double hz = max[2] - min_z;

                int n_elements = nx * ny * nz;
                int n_points = (nx + 1) * (ny + 1) * (nz + 1);

                mesh.set_n_parts(1 + 6);
                mesh.set_point_dims(n_points, 3);
                mesh.parts[0].set_elem_dims(n_elements, 8);
                mesh.parts[0].elem_type = 8;
                mesh.parts[0].part_id = 0;

                for (int i = 1; i <= 6; ++i) {
                    mesh.parts[i].elem_type = 4;
                    mesh.parts[i].part_id = i;
                }

                mesh.parts[0].name = "bulk";
                mesh.parts[1].name = "left";
                mesh.parts[2].name = "bottom";
                mesh.parts[3].name = "back";
                mesh.parts[4].name = "right";
                mesh.parts[5].name = "top";
                mesh.parts[6].name = "front";

                // Left
                mesh.parts[1].set_elem_dims(ny * nz, 4);
                //  Right
                mesh.parts[4].set_elem_dims(ny * nz, 4);

                // Bottom
                mesh.parts[2].set_elem_dims(nx * nz, 4);

                // Top
                mesh.parts[5].set_elem_dims(nx * nz, 4);

                // Back
                mesh.parts[3].set_elem_dims(nx * ny, 4);

                // Front
                mesh.parts[6].set_elem_dims(nx * ny, 4);

                auto points = mesh.points;

                Kokkos::parallel_for(
                    "Cube::build::BoundaryFlags",
                    Range3({0, 0, 0}, {nx + 1, ny + 1, nz + 1}),
                    KOKKOS_LAMBDA(int i, int j, int k) {
                        int inode = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;

                        points(inode, 0) = min_x + (double)i * hx;
                        points(inode, 1) = min_y + (double)j * hy;
                        points(inode, 2) = min_z + (double)k * hz;
                    });

                auto elem_to_node = mesh.parts[0].elem_to_node;

                Kokkos::parallel_for(
                    "Cube::build::elements", Range3({0, 0, 0}, {nx, ny, nz}), KOKKOS_LAMBDA(int i, int j, int k) {
                        int ielem = k * ny * nx + j * nx + i;
                        elem_to_node(ielem, 0) = (ny + 1) * (nx + 1) * k + (nx + 1) * j + i;
                        elem_to_node(ielem, 1) = (ny + 1) * (nx + 1) * k + (nx + 1) * j + i + 1;
                        elem_to_node(ielem, 2) = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + i + 1;
                        elem_to_node(ielem, 3) = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + i;
                        elem_to_node(ielem, 4) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * j + i;
                        elem_to_node(ielem, 5) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * j + i + 1;
                        elem_to_node(ielem, 6) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * (j + 1) + i + 1;
                        elem_to_node(ielem, 7) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * (j + 1) + i;
                    });

                auto left = mesh.parts[1].elem_to_node;
                auto bottom = mesh.parts[2].elem_to_node;
                auto back = mesh.parts[3].elem_to_node;
                auto right = mesh.parts[4].elem_to_node;
                auto top = mesh.parts[5].elem_to_node;
                auto front = mesh.parts[6].elem_to_node;

                Kokkos::parallel_for(
                    "Cube::build::left", Range2({0, 0}, {ny, nz}), KOKKOS_LAMBDA(int j, int k) {
                        int ielem = k * ny + j;
                        left(ielem, 0) = (ny + 1) * (nx + 1) * k + (nx + 1) * j + 0;
                        left(ielem, 1) = (ny + 1) * (nx + 1) * k + (nx + 1) * j + 0 + 1;
                        left(ielem, 2) = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + 0 + 1;
                        left(ielem, 3) = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + 0;
                    });

                Kokkos::parallel_for(
                    "Cube::build::right", Range2({0, 0}, {ny, nz}), KOKKOS_LAMBDA(int j, int k) {
                        int ielem = k * ny + j;
                        right(ielem, 0) = (ny + 1) * (nx + 1) * k + (nx + 1) * j + nx;
                        right(ielem, 1) = (ny + 1) * (nx + 1) * k + (nx + 1) * j + nx + 1;
                        right(ielem, 2) = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + nx + 1;
                        right(ielem, 3) = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + nx;
                    });

                Kokkos::parallel_for(
                    "Cube::build::back", Range2({0, 0}, {nx, ny}), KOKKOS_LAMBDA(int i, int j) {
                        int ielem = j * nx + i;
                        bottom(ielem, 0) = (ny + 1) * (nx + 1) * 0 + (nx + 1) * j + i;
                        bottom(ielem, 1) = (ny + 1) * (nx + 1) * 0 + (nx + 1) * j + i + 1;
                        bottom(ielem, 2) = (ny + 1) * (nx + 1) * 0 + (nx + 1) * (j + 1) + i + 1;
                        bottom(ielem, 3) = (ny + 1) * (nx + 1) * 0 + (nx + 1) * (j + 1) + i;
                    });

                Kokkos::parallel_for(
                    "Cube::build::front", Range2({0, 0}, {nx, ny}), KOKKOS_LAMBDA(int i, int j) {
                        int ielem = j * nx + i;
                        top(ielem, 0) = (ny + 1) * (nx + 1) * nz + (nx + 1) * j + i;
                        top(ielem, 1) = (ny + 1) * (nx + 1) * nz + (nx + 1) * j + i + 1;
                        top(ielem, 2) = (ny + 1) * (nx + 1) * nz + (nx + 1) * (j + 1) + i + 1;
                        top(ielem, 3) = (ny + 1) * (nx + 1) * nz + (nx + 1) * (j + 1) + i;
                    });

                Kokkos::parallel_for(
                    "Cube::build::bottom", Range2({0, 0}, {nx, nz}), KOKKOS_LAMBDA(int i, int k) {
                        int ielem = k * nx + i;
                        top(ielem, 0) = (ny + 1) * (nx + 1) * nz + (nx + 1) * 0 + i;
                        top(ielem, 1) = (ny + 1) * (nx + 1) * nz + (nx + 1) * 0 + i + 1;
                        top(ielem, 2) = (ny + 1) * (nx + 1) * nz + (nx + 1) * (0 + 1) + i + 1;
                        top(ielem, 3) = (ny + 1) * (nx + 1) * nz + (nx + 1) * (0 + 1) + i;
                    });

                Kokkos::parallel_for(
                    "Cube::build::top", Range2({0, 0}, {nx, nz}), KOKKOS_LAMBDA(int i, int k) {
                        int ielem = k * nx + i;
                        top(ielem, 0) = (ny + 1) * (nx + 1) * nz + (nx + 1) * ny + i;
                        top(ielem, 1) = (ny + 1) * (nx + 1) * nz + (nx + 1) * ny + i + 1;
                        top(ielem, 2) = (ny + 1) * (nx + 1) * nz + (nx + 1) * (ny + 1) + i + 1;
                        top(ielem, 3) = (ny + 1) * (nx + 1) * nz + (nx + 1) * (ny + 1) + i;
                    });

                c.stop();

                if (debug) {
                    out() << "Cube::build(...)\nMesh with " << n_points << " points\n";
                    out() << "parts:\n";
                    for (auto &p : mesh.parts) {
                        out() << "\t" << p.part_id << ") " << p.elem_to_node.extent(0) << "\t"
                              << (p.elem_type == 8 ? "HEX8 " : "QUAD4") << " elements\t (" << p.name << ")\n";
                    }

                    out() << c << "\n";
                }
            }

            int n[3] = {1, 1, 1};
            double min[3] = {0., 0., 0.};
            double max[3] = {1., 1., 1.};
            bool debug{false};
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_CUBE_HPP

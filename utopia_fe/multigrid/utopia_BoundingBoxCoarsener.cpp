#include "utopia_BoundingBoxCoarsener.hpp"

// #include "moonolith_profiler.hpp"
// #include "moonolith_redistribute.hpp"
// #include "moonolith_tree.hpp"
// #include "moonolith_n_tree_mutator_factory.hpp"
// #include "moonolith_n_tree_with_span_mutator_factory.hpp"
// #include "moonolith_n_tree_with_tagsmutator_factory.hpp"
// #include "moonolith_sparse_matrix.hpp"
// #include "moonolith_communicator.hpp"

#include "moonolith_aabb.hpp"
#include "par_moonolith.hpp"

#include "utopia_Contact.hpp"
#include "utopia_Socket.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_libmesh.hpp"

#include "libmesh/explicit_system.h"
#include "libmesh/libmesh_version.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/serial_mesh.h"

#include <array>
#include <map>

namespace utopia {
    class BoundingBoxCoarsener::Impl {
    public:
        template <moonolith::Integer Dimension>
        static void expand_box(const libMesh::MeshBase &mesh,
                               const libMesh::Elem &e,
                               moonolith::AABB<Dimension, double> &box) {
            std::array<double, Dimension> p_a;
            for (libMesh::dof_id_type i = 0; i < e.n_nodes(); ++i) {
                const libMesh::Point &p = mesh.node(e.node_id(i));
                for (int d = 0; d < Dimension; ++d) {
                    p_a[d] = p(d);
                }

                box += p_a;
            }
        }

        template <moonolith::Integer Dimension>
        static void synch_box(moonolith::Communicator &comm, moonolith::AABB<Dimension, double> &box) {
            std::array<double, 2 * Dimension> min_max;

            for (moonolith::Integer d = 0; d < Dimension; ++d) {
                min_max[d] = box.min(d);
                min_max[Dimension + d] = -box.max(d);
            }

            comm.all_reduce(&min_max[0], min_max.size(), moonolith::MPIMin());

            for (moonolith::Integer d = 0; d < Dimension; ++d) {
                box.set_min(d, min_max[d]);
                box.set_max(d, -min_max[Dimension + d]);
            }
        }

        static void synch_set(moonolith::Communicator &comm, std::set<int> &s) {
            moonolith::ByteOutputBuffer send_buffer;
            std::vector<moonolith::ByteInputBuffer> recv_buffer;

            moonolith::write_set(s, send_buffer);
            comm.unstructured_all_gather(send_buffer, recv_buffer, true);

            for (auto &rb : recv_buffer) {
                if (rb.good()) {
                    moonolith::read_set(s, rb);
                }
            }
        }

        static void synch_add_map(moonolith::Communicator &comm, std::map<int, long> &m) {
            moonolith::ByteOutputBuffer send_buffer;
            std::vector<moonolith::ByteInputBuffer> recv_buffer;

            moonolith::write_map(m, send_buffer);
            comm.unstructured_all_gather(send_buffer, recv_buffer, true);

            std::map<int, long> buffer;
            for (auto &rb : recv_buffer) {
                if (rb.good()) {
                    buffer.clear();
                    moonolith::read_map(buffer, rb);

                    for (auto &b : buffer) {
                        m[b.first] += b.second;
                    }
                }
            }
        }

        inline bool init_2d(const int n_coarsening_levels, const libMesh::MeshBase &mesh) {
            return init_aux<2>(n_coarsening_levels, mesh);
        }

        inline bool init_3d(const int n_coarsening_levels, const libMesh::MeshBase &mesh) {
            return init_aux<3>(n_coarsening_levels, mesh);
        }

        static void set_subdomain_ids(const int tag, libMesh::MeshBase &mesh) {
            for (auto it = elements_begin(mesh); it != elements_end(mesh); ++it) {
                (*it)->subdomain_id() = tag;
            }
        }

        static void build_mesh_from_box(const moonolith::AABB<2, double> &box,
                                        const moonolith::Integer n_nodes,
                                        const int n_coarsening_levels,
                                        libMesh::UnstructuredMesh &mesh) {
            std::array<double, 2> r{box.max(0) - box.min(0), box.max(1) - box.min(1)};

            const int n_segments = std::max(2, int(std::round(std::sqrt(n_nodes / std::pow(4, n_coarsening_levels)))));
            const double max_r = std::max(r[0], r[1]);

            const double aspect_ratio_x = r[0] / max_r;
            const double aspect_ratio_y = r[1] / max_r;

            const int nx = std::max(int(std::round(n_segments * aspect_ratio_x)), 2);
            const int ny = std::max(int(std::round(n_segments * aspect_ratio_y)), 2);

            libMesh::MeshTools::Generation::build_square(
                mesh, nx, ny, box.min(0), box.max(0), box.min(1), box.max(1), libMesh::QUAD4);
        }

        static void build_mesh_from_box(const moonolith::AABB<3, double> &box,
                                        const moonolith::Integer n_elem,
                                        const int n_coarsening_levels,
                                        libMesh::UnstructuredMesh &mesh) {
            std::array<double, 3> r{box.max(0) - box.min(0), box.max(1) - box.min(1), box.max(2) - box.min(2)};

            const int n_segments = std::max(2, int(std::round(std::cbrt(n_elem / std::pow(8, n_coarsening_levels)))));
            const double max_r = std::max(std::max(r[0], r[1]), r[2]);

            const double aspect_ratio_x = r[0] / max_r;
            const double aspect_ratio_y = r[1] / max_r;
            const double aspect_ratio_z = r[2] / max_r;

            const int nx = std::max(int(std::round(n_segments * aspect_ratio_x)), 2);
            const int ny = std::max(int(std::round(n_segments * aspect_ratio_y)), 2);
            const int nz = std::max(int(std::round(n_segments * aspect_ratio_z)), 2);

            libMesh::MeshTools::Generation::build_cube(mesh,
                                                       nx,
                                                       ny,
                                                       nz,
                                                       box.min(0),
                                                       box.max(0),
                                                       box.min(1),
                                                       box.max(1),
                                                       box.min(2),
                                                       box.max(2),
                                                       libMesh::HEX8);
        }

        static void concatenate_distributed_meshes(std::vector<std::shared_ptr<libMesh::UnstructuredMesh>> &sub_meshes,
                                                   const int tag_id_offset,
                                                   libMesh::MeshBase &coarse_mesh) {
            coarse_mesh.clear();
            coarse_mesh.set_mesh_dimension(sub_meshes.front()->mesh_dimension());
            coarse_mesh.set_spatial_dimension(sub_meshes.front()->mesh_dimension());

            std::vector<long> node_offsets(sub_meshes.size() + 1, 0);
            std::vector<long> element_offsets(sub_meshes.size() + 1, 0);

            {
                std::size_t i = 0;
                for (auto &m_ptr : sub_meshes) {
                    node_offsets[i + 1] = m_ptr->n_nodes() + node_offsets[i];
                    element_offsets[i + 1] = m_ptr->n_active_elem() + element_offsets[i];
                    ++i;
                }
            }

            moonolith::Communicator comm(coarse_mesh.comm().get());
            for (int i = 0; i < comm.size(); ++i) {
                comm.barrier();
                if (i != comm.rank()) continue;

                std::size_t sub_index = 0;
                for (auto m_ptr : sub_meshes) {
                    std::cout << moonolith::Communicator() << "\n------------------------------" << std::endl;

                    // for(auto n_it = m_ptr->local_nodes_begin(); n_it != m_ptr->local_nodes_end(); ++n_it) {
                    for (auto n_it = m_ptr->active_nodes_begin(); n_it != m_ptr->active_nodes_end(); ++n_it) {
                        auto n_ptr = coarse_mesh.add_point(
                            **n_it, node_offsets[sub_index] + (*n_it)->id(), (*n_it)->processor_id());
                        n_ptr->set_unique_id() = (node_offsets[sub_index] + (*n_it)->unique_id());

                        std::cout << "node : " << (*n_it)->id() << " -> " << n_ptr->id() << std::endl;
                    }

                    std::cout << "-----------------------------" << std::endl;

                    for (auto e_it = m_ptr->active_local_elements_begin(); e_it != m_ptr->active_local_elements_end();
                         ++e_it) {
                        auto &e = **e_it;
                        auto elem = libMesh::Elem::build(libMesh::ElemType(e.type())).release();

                        for (int ii = 0; ii != e.n_nodes(); ++ii) {
                            elem->set_node(ii) = coarse_mesh.node_ptr(e.node_id(ii) + node_offsets[sub_index]);
                        }

                        // elem->set_id(e.id() + element_offsets[sub_index]);
                        elem->subdomain_id() = e.subdomain_id() + tag_id_offset;
                        elem->processor_id() = e.processor_id();
                        elem->set_unique_id() = element_offsets[sub_index] + e.unique_id();

                        // std::cout << moonolith::Communicator() << " ";

                        coarse_mesh.add_elem(elem);

                        std::cout << e.unique_id() << " " << elem->id() << "( " << e.id()
                                  << " ) : " << elem->subdomain_id() << "->" << (e.subdomain_id()) << std::endl;
                    }

                    // node_id_offset += node_count;
                    ++sub_index;

                    std::cout << "-----------------------------" << std::endl;
                }
            }

            comm.barrier();

            coarse_mesh.prepare_for_use(/*skip_renumber =*/false);
        }

        static void concatenate_meshes(std::vector<std::shared_ptr<libMesh::UnstructuredMesh>> &sub_meshes,
                                       const int tag_id_offset,
                                       libMesh::MeshBase &coarse_mesh) {
            coarse_mesh.clear();
            coarse_mesh.set_mesh_dimension(sub_meshes.front()->mesh_dimension());
            coarse_mesh.set_spatial_dimension(sub_meshes.front()->mesh_dimension());

            moonolith::Communicator comm(coarse_mesh.comm().get());

            for (int i = 0; i < comm.size(); ++i) {
                comm.barrier();
                if (i != comm.rank()) continue;

                auto node_id_offset = 0;
                std::size_t element_id = 0;
                for (auto m_ptr : sub_meshes) {
                    std::cout << moonolith::Communicator() << "\n------------------------------" << std::endl;

                    std::size_t node_count = 0;
                    for (auto n_it = m_ptr->active_nodes_begin(); n_it != m_ptr->active_nodes_end(); ++n_it) {
                        auto n_ptr =
                            coarse_mesh.add_point(**n_it, node_id_offset + (*n_it)->id(), coarse_mesh.processor_id());
                        ++node_count;
                    }

                    for (auto e_it = m_ptr->active_elements_begin(); e_it != m_ptr->active_elements_end(); ++e_it) {
                        auto &e = **e_it;
                        auto elem = libMesh::Elem::build(libMesh::ElemType(e.type())).release();

                        for (int ii = 0; ii != e.n_nodes(); ++ii) {
                            elem->set_node(ii) = coarse_mesh.node_ptr(e.node_id(ii) + node_id_offset);
                        }

                        elem->set_id(element_id++);
                        elem->subdomain_id() = e.subdomain_id() + tag_id_offset;

                        // std::cout << moonolith::Communicator() << " ";
                        std::cout << e.unique_id() << " " << elem->id() << "( " << e.id()
                                  << " ) : " << elem->subdomain_id() << "->" << (e.subdomain_id()) << std::endl;
                        coarse_mesh.add_elem(elem);
                    }

                    node_id_offset += node_count;
                }

                std::cout << "\n------------------------------" << std::endl;
            }

            coarse_mesh.prepare_for_use(/*skip_renumber =*/false);
        }

        template <moonolith::Integer Dimension>
        bool init_aux(const int n_coarsening_levels, const libMesh::MeshBase &mesh) {
            std::map<int, moonolith::AABB<Dimension, double>> boxes;
            std::map<int, long> n_elems;

            moonolith::Communicator comm(mesh.comm().get());

            for (auto it = elements_begin(mesh); it != elements_end(mesh); ++it) {
                auto &e = **it;
                auto &box = boxes[e.subdomain_id()];
                expand_box(mesh, e, box);
                ++n_elems[e.subdomain_id()];
            }

            std::set<int> box_ids;
            for (auto &b : boxes) {
                box_ids.insert(b.first);
            }

            synch_set(comm, box_ids);
            synch_add_map(comm, n_elems);

            if (box_ids.empty()) {
                assert(false);
                std::cerr << "something went wrong" << std::endl;
                return false;
            }

            const int tag_id_offset = *box_ids.rbegin() + 1;

            tags.reserve(box_ids.size());

            for (auto id : box_ids) {
                tags.emplace_back(tag_id_offset + id, id);
            }

            for (auto id : box_ids) {
                // allocate additional boxes
                boxes[id];
            }

            for (auto &b : boxes) {
                synch_box(comm, b.second);
            }

            std::vector<std::shared_ptr<libMesh::UnstructuredMesh>> sub_meshes;

            for (auto &b : boxes) {
                auto m_ptr = std::make_shared<libMesh::SerialMesh>(mesh.comm());
                sub_meshes.push_back(m_ptr);
                build_mesh_from_box(b.second, n_elems[b.first], n_coarsening_levels, *m_ptr);
                set_subdomain_ids(b.first, *m_ptr);
            }

            if (sub_meshes.size() == 1) {
                coarse_mesh = sub_meshes.front();
            } else {
                coarse_mesh = std::make_shared<libMesh::DistributedMesh>(mesh.comm());
                concatenate_distributed_meshes(sub_meshes, tag_id_offset, *coarse_mesh);
            }

            return true;
        }

        std::shared_ptr<libMesh::UnstructuredMesh> coarse_mesh;
        std::vector<std::pair<int, int>> tags;
    };

    BoundingBoxCoarsener::~BoundingBoxCoarsener() {}

    BoundingBoxCoarsener::BoundingBoxCoarsener() : impl_(std::make_shared<BoundingBoxCoarsener::Impl>()) {}

    bool BoundingBoxCoarsener::init(const int n_coarsening_levels, const libMesh::MeshBase &mesh) {
        switch (mesh.mesh_dimension()) {
            case 2: {
                return impl_->init_2d(n_coarsening_levels, mesh);
            }

            case 3: {
                return impl_->init_3d(n_coarsening_levels, mesh);
            }

            default: { return false; }
        }
    }

    void BoundingBoxCoarsener::describe(std::ostream &os) const {
        if (impl_->coarse_mesh) {
            impl_->coarse_mesh->print_info(os);

            if (mpi_world_size() == 1) {
                plot_mesh(*impl_->coarse_mesh, "mesh/coarse");
            }
        }
    }

    const std::vector<std::pair<int, int>> &BoundingBoxCoarsener::get_tags() const { return impl_->tags; }

    const std::shared_ptr<libMesh::UnstructuredMesh> &BoundingBoxCoarsener::get_mesh() const {
        return impl_->coarse_mesh;
    }
}  // namespace utopia

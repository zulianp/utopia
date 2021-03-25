#include "utopia_moonolith_stk.hpp"

// Utopia fe
#include "utopia_moonolith_Mesh.hpp"
#include "utopia_stk_Mesh.hpp"

#include "utopia_moonolith_FunctionSpace.hpp"
#include "utopia_stk_Commons.hpp"
#include "utopia_stk_FunctionSpace.hpp"

// Moonolith
#include "moonolith_elem_type.hpp"
#include "moonolith_mesh.hpp"

// Stk
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>

// Std lib
#include <cassert>

namespace utopia {

    template <int Dim>
    class ConvertMeshSTK2Moonolith {
    public:
        using Bucket_t = ::stk::mesh::Bucket;
        using BucketVector_t = ::stk::mesh::BucketVector;
        using Entity_t = ::stk::mesh::Entity;
        using Scalar_t = Traits<utopia::stk::Mesh>::Scalar;
        using Size_t = Traits<utopia::stk::Mesh>::SizeType;
        using MoonolithMesh_t = ::moonolith::Mesh<Scalar_t, Dim>;
        using MoonolithFunctionSpace_t = ::moonolith::FunctionSpace<::moonolith::Mesh<Scalar_t, Dim>>;
        using MetaData_t = ::stk::mesh::MetaData;
        using BulkData_t = ::stk::mesh::BulkData;

        static int convert_to_manifold_dim(::moonolith::ElemType type) {
            switch (type) {
                case ::moonolith::NODE1:
                    return 0;
                case ::moonolith::EDGE2:
                    return 1;
                case ::moonolith::TRI3:
                case ::moonolith::QUAD4:
                    return 2;
                case ::moonolith::HEX8:
                case ::moonolith::HEX27:
                case ::moonolith::TET4:
                    return 3;
                default: {
                    assert(false);
                    return -1;
                }
            }
        }

        static ::moonolith::ElemType convert_elem_type(::stk::topology::topology_t topo) {
            switch (topo) {
                case ::stk::topology::NODE:
                    return ::moonolith::NODE1;
                case ::stk::topology::LINE_2:
                    return ::moonolith::EDGE2;
                case ::stk::topology::TRI_3:
                case ::stk::topology::TRI_3_2D:
                case ::stk::topology::SHELL_TRI_3:
                    return ::moonolith::TRI3;
                case ::stk::topology::QUAD_4:
                case ::stk::topology::QUAD_4_2D:
                case ::stk::topology::SHELL_QUAD_4:
                    return ::moonolith::QUAD4;
                case ::stk::topology::HEX_8:
                    return ::moonolith::HEX8;
                case ::stk::topology::HEX_27:
                    return ::moonolith::HEX27;
                case ::stk::topology::TET_4:
                    return ::moonolith::TET4;
                default: {
                    assert(false);
                    return ::moonolith::INVALID;
                }
            }
        }

        static void copy_meta_info(const utopia::stk::FunctionSpace &in, utopia::moonolith::FunctionSpace &out) {
            auto m_space = out.raw_type<Dim>();
            m_space->dof_map().set_n_local_dofs(in.mesh().n_local_nodes());
            m_space->dof_map().set_n_dofs(in.mesh().n_nodes());
        }

        static void apply(const utopia::stk::Mesh &in, utopia::moonolith::Mesh &out) {
            assert(in.comm().size() == 1 && "IMPLEMENT PARALLEL VERSION");

            auto &meta_data = in.meta_data();
            auto &bulk_data = in.bulk_data();

            // const Size_t n_elements = in.n_elements();
            // const Size_t n_nodes = in.n_nodes();
            const Size_t n_local_elements = in.n_local_elements();
            const Size_t n_local_nodes = in.n_local_nodes();

            auto m_mesh = std::make_shared<MoonolithMesh_t>(in.comm().raw_comm());

            m_mesh->resize(n_local_elements, n_local_nodes);

            ::stk::mesh::Selector s_universal = meta_data.universal_part();
            const BucketVector_t &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, s_universal);

            ::stk::mesh::FieldBase *coords = ::stk::mesh::get_field_by_name("coordinates", meta_data);

            for (const auto &ib : node_buckets) {
                const Bucket_t &b = *ib;
                const Bucket_t::size_type length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    Entity_t node = b[k];
                    auto moonolith_index = utopia::stk::convert_entity_to_index(node);
                    auto &p = m_mesh->node(moonolith_index);

                    const Scalar_t *points = (const Scalar_t *)::stk::mesh::field_data(*coords, node);

                    for (int d = 0; d < Dim; ++d) {
                        p[d] = points[d];
                    }
                }
            }

            const BucketVector_t &elem_buckets = bulk_data.get_buckets(::stk::topology::ELEMENT_RANK, s_universal);

            int manifold_dim = -1;

            for (const auto &ib : elem_buckets) {
                const Bucket_t &b = *ib;
                const Bucket_t::size_type length = b.size();
                auto moonolith_type = convert_elem_type(b.topology());
                manifold_dim = std::max(manifold_dim, convert_to_manifold_dim(moonolith_type));

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    // get the current node entity and extract the id to fill it into the field
                    Entity_t elem = b[k];
                    const Size_t elem_idx = utopia::stk::convert_entity_to_index(elem) - n_local_nodes;
                    const Size_t n_nodes = bulk_data.num_nodes(elem);

                    auto &e = m_mesh->elem(elem_idx);
                    e.type = moonolith_type;
                    e.block = 1;              // elem.subdomain_id();
                    e.is_affine = true;       /// elem.has_affine_map();
                    e.global_idx = elem_idx;  // elem.unique_id();
                    e.nodes.resize(n_nodes);

                    auto node_ids = bulk_data.begin_nodes(elem);

                    for (Size_t i = 0; i < n_nodes; ++i) {
                        e.nodes[i] = utopia::stk::convert_entity_to_index(node_ids[i]);
                    }
                }
            }

            m_mesh->set_manifold_dim(manifold_dim);
            m_mesh->finalize();
            out.wrap(m_mesh);
        }

        static void extract_surface(const utopia::stk::Mesh &in, utopia::moonolith::Mesh &out) {
            auto &meta_data = in.meta_data();
            auto &bulk_data = in.bulk_data();

            ::stk::mesh::Selector s_universal = meta_data.universal_part();

            ::stk::mesh::FieldBase *coords = ::stk::mesh::get_field_by_name("coordinates", meta_data);

            auto topo = ::stk::topology::FACE_RANK;

            if (in.spatial_dimension() == 2) {
                topo = ::stk::topology::EDGE_RANK;
            } else {
                assert(in.spatial_dimension() == 3);
            }

            const BucketVector_t &elem_buckets = bulk_data.get_buckets(topo, s_universal);

            Bucket_t::size_type n_selected_elements = 0;
            Bucket_t::size_type n_selected_nodes = 0;
            Bucket_t::size_type n_local_nodes = in.n_local_nodes();

            std::vector<SizeType> node_mapping(n_local_nodes, -1);

            for (const auto &b_ptr : elem_buckets) {
                auto &b = *b_ptr;

                const Bucket_t::size_type length = b.size();

                n_selected_elements += length;

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    Entity_t elem = b[k];

                    const Size_t n_nodes = bulk_data.num_nodes(elem);
                    auto node_ids = bulk_data.begin_nodes(elem);

                    for (Size_t i = 0; i < n_nodes; ++i) {
                        auto idx = utopia::stk::convert_entity_to_index(node_ids[i]);

                        if (node_mapping[idx] == -1) {
                            node_mapping[idx] = n_selected_nodes++;
                        }
                    }
                }
            }

            // std::cout << "Num elem: " << n_selected_elements << "\n";

            ////////////////////////////////////////////////////////////////
            auto m_mesh = std::make_shared<MoonolithMesh_t>(in.comm().raw_comm());

            m_mesh->resize(n_selected_elements, n_selected_nodes);

            ////////////////////////////////////////////////////////////////
            // Nodes
            ////////////////////////////////////////////////////////////////

            const BucketVector_t &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, s_universal);

            for (const auto &ib : node_buckets) {
                const Bucket_t &b = *ib;
                const Bucket_t::size_type length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    Entity_t node = b[k];
                    auto moonolith_index = node_mapping[utopia::stk::convert_entity_to_index(node)];

                    if (moonolith_index == -1) continue;

                    assert(moonolith_index < n_selected_nodes);

                    auto &p = m_mesh->node(moonolith_index);

                    const Scalar_t *points = (const Scalar_t *)::stk::mesh::field_data(*coords, node);

                    for (int d = 0; d < Dim; ++d) {
                        p[d] = points[d];
                    }
                }
            }

            ////////////////////////////////////////////////////////////////
            // Elements
            ////////////////////////////////////////////////////////////////

            SizeType selected_elem_idx = 0;
            for (const auto &ib : elem_buckets) {
                const Bucket_t &b = *ib;

                auto moonolith_type = convert_elem_type(b.topology());

                const Bucket_t::size_type length = b.size();

                // std::cout << "Bucket with " << length << " ";
                // std::cout << b.topology() << " -> ";
                // std::cout << moonolith_type << "\n";

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    Entity_t elem = b[k];
                    const Size_t elem_idx = utopia::stk::convert_entity_to_index(elem) - n_local_nodes;
                    const Size_t n_nodes = bulk_data.num_nodes(elem);

                    auto &e = m_mesh->elem(selected_elem_idx++);
                    e.type = moonolith_type;
                    e.block = 1;              // elem.subdomain_id();
                    e.is_affine = true;       /// elem.has_affine_map();
                    e.global_idx = elem_idx;  // elem.unique_id();
                    e.nodes.resize(n_nodes);

                    auto node_ids = bulk_data.begin_nodes(elem);

                    for (Size_t i = 0; i < n_nodes; ++i) {
                        e.nodes[i] = node_mapping[utopia::stk::convert_entity_to_index(node_ids[i])];
                    }
                }
            }

            m_mesh->set_manifold_dim(in.spatial_dimension() - 1);
            m_mesh->finalize();

            out.wrap(m_mesh);
        }

        static void extract_trace_space(const utopia::stk::FunctionSpace &in, utopia::moonolith::FunctionSpace &out) {
            auto m_mesh = std::make_shared<utopia::moonolith::Mesh>(in.comm());
            extract_surface(in.mesh(), *m_mesh);
            out.init(m_mesh, false);

            auto &out_space = *out.raw_type<Dim>();

            auto topo = ::stk::topology::FACE_RANK;

            if (in.mesh().spatial_dimension() == 2) {
                topo = ::stk::topology::EDGE_RANK;
            }

            auto &meta_data = in.mesh().meta_data();
            auto &bulk_data = in.mesh().bulk_data();

            ::stk::mesh::Selector s_universal = meta_data.universal_part();
            const BucketVector_t &elem_buckets = bulk_data.get_buckets(topo, s_universal);

            // Bucket_t::size_type n_selected_elements = 0;
            // Bucket_t::size_type n_selected_nodes = 0;
            Bucket_t::size_type n_local_nodes = in.mesh().n_local_nodes();
            Bucket_t::size_type n_local_dofs = in.n_local_dofs();
            Bucket_t::size_type n_dofs = in.n_dofs();

            auto &out_dof_map = out_space.dof_map();
            out_dof_map.resize(m_mesh->n_local_elements());
            out_dof_map.set_n_local_dofs(n_local_dofs);
            out_dof_map.set_n_dofs(n_dofs);
            // out_dof_map.set_max_nnz(max_nnz_x_row(dof_map));

            SizeType selected_elem_idx = 0;
            for (const auto &ib : elem_buckets) {
                const Bucket_t &b = *ib;

                auto moonolith_type = convert_elem_type(b.topology());

                const Bucket_t::size_type length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    Entity_t elem = b[k];
                    const Size_t elem_idx = utopia::stk::convert_entity_to_index(elem) - n_local_nodes;
                    const Size_t n_nodes = bulk_data.num_nodes(elem);

                    auto &dof_object = out_dof_map.dof_object(selected_elem_idx++);
                    dof_object.type = moonolith_type;
                    dof_object.block = 1;              // elem.subdomain_id();
                    dof_object.global_idx = elem_idx;  // elem.unique_id();
                    dof_object.dofs.resize(n_nodes);
                    dof_object.element_dof = utopia::stk::convert_entity_to_index(elem);

                    auto node_ids = bulk_data.begin_nodes(elem);

                    for (Size_t i = 0; i < n_nodes; ++i) {
                        dof_object.dofs[i] = utopia::stk::convert_entity_to_index(node_ids[i]);
                    }
                }
            }

            assert(out_space.dof_map().is_valid());
        }
    };

    void ConvertMesh<utopia::stk::Mesh, utopia::moonolith::Mesh>::apply(const utopia::stk::Mesh &in,
                                                                        utopia::moonolith::Mesh &out) {
        auto &meta_data = in.meta_data();
        // auto &bulk_data = in.bulk_data();

        const int dim = meta_data.spatial_dimension();

        switch (dim) {
            case 1: {
                ConvertMeshSTK2Moonolith<1>::apply(in, out);
                break;
            }

            case 2: {
                ConvertMeshSTK2Moonolith<2>::apply(in, out);
                break;
            }

            case 3: {
                ConvertMeshSTK2Moonolith<3>::apply(in, out);
                break;
            }
            default: {
                Utopia::Abort();
            }
        }
    }

    void ConvertFunctionSpace<utopia::stk::FunctionSpace, utopia::moonolith::FunctionSpace>::apply(
        const utopia::stk::FunctionSpace &in,
        utopia::moonolith::FunctionSpace &out) {
        auto m_mesh = std::make_shared<utopia::moonolith::Mesh>(in.comm());
        convert_mesh(in.mesh(), *m_mesh);

        out.init(m_mesh);

        const int dim = in.mesh().spatial_dimension();

        switch (dim) {
            case 1: {
                ConvertMeshSTK2Moonolith<1>::copy_meta_info(in, out);
                break;
            }

            case 2: {
                ConvertMeshSTK2Moonolith<2>::copy_meta_info(in, out);
                break;
            }

            case 3: {
                ConvertMeshSTK2Moonolith<3>::copy_meta_info(in, out);
                break;
            }
            default: {
                Utopia::Abort();
            }
        }
    }

    void ExtractSurface<utopia::stk::Mesh, utopia::moonolith::Mesh>::apply(const utopia::stk::Mesh &in,
                                                                           utopia::moonolith::Mesh &out) {
        auto &meta_data = in.meta_data();

        const int dim = meta_data.spatial_dimension();

        switch (dim) {
            case 1: {
                ConvertMeshSTK2Moonolith<1>::extract_surface(in, out);
                break;
            }

            case 2: {
                ConvertMeshSTK2Moonolith<2>::extract_surface(in, out);
                break;
            }

            case 3: {
                ConvertMeshSTK2Moonolith<3>::extract_surface(in, out);
                break;
            }
            default: {
                Utopia::Abort();
            }
        }
    }

    void ExtractTraceSpace<utopia::stk::FunctionSpace, utopia::moonolith::FunctionSpace>::apply(
        const utopia::stk::FunctionSpace &in,
        utopia::moonolith::FunctionSpace &out) {
        auto &meta_data = in.mesh().meta_data();

        const int dim = meta_data.spatial_dimension();

        switch (dim) {
            case 1: {
                ConvertMeshSTK2Moonolith<1>::extract_trace_space(in, out);
                break;
            }

            case 2: {
                ConvertMeshSTK2Moonolith<2>::extract_trace_space(in, out);
                break;
            }

            case 3: {
                ConvertMeshSTK2Moonolith<3>::extract_trace_space(in, out);
                break;
            }
            default: {
                Utopia::Abort();
            }
        }
    }

}  // namespace utopia

// enum topology_t
// {
//     INVALID_TOPOLOGY
//   , BEGIN_TOPOLOGY
//   //NODE_RANK
//   , NODE = BEGIN_TOPOLOGY
//   //EDGE_RANK
//   , LINE_2
//   , LINE_3
//   //FACE_RANK
//   , TRI_3, TRIANGLE_3 = TRI_3
//   , TRI_4, TRIANGLE_4 = TRI_4
//   , TRI_6, TRIANGLE_6 = TRI_6
//   , QUAD_4, QUADRILATERAL_4 = QUAD_4
//   , QUAD_6, QUADRILATERAL_6 = QUAD_6
//   , QUAD_8, QUADRILATERAL_8 = QUAD_8
//   , QUAD_9, QUADRILATERAL_9 = QUAD_9
//   //ELEMENT_RANK
//   , PARTICLE, BEGIN_ELEMENT_RANK = PARTICLE
//   , LINE_2_1D
//   , LINE_3_1D
//   , BEAM_2
//   , BEAM_3
//   , SHELL_LINE_2
//   , SHELL_LINE_3
//   , SPRING_2
//   , SPRING_3
//   , TRI_3_2D, TRIANGLE_3_2D = TRI_3_2D
//   , TRI_4_2D, TRIANGLE_4_2D = TRI_4_2D
//   , TRI_6_2D, TRIANGLE_6_2D = TRI_6_2D
//   , QUAD_4_2D, QUADRILATERAL_4_2D = QUAD_4_2D
//   , QUAD_8_2D, QUADRILATERAL_8_2D = QUAD_8_2D
//   , QUAD_9_2D, QUADRILATERAL_9_2D = QUAD_9_2D
//   , SHELL_TRI_3, SHELL_TRIANGLE_3 = SHELL_TRI_3
//   , SHELL_TRI_4, SHELL_TRIANGLE_4 = SHELL_TRI_4
//   , SHELL_TRI_6, SHELL_TRIANGLE_6 = SHELL_TRI_6
//   , SHELL_QUAD_4, SHELL_QUADRILATERAL_4 = SHELL_QUAD_4
//   , SHELL_QUAD_8, SHELL_QUADRILATERAL_8 = SHELL_QUAD_8
//   , SHELL_QUAD_9, SHELL_QUADRILATERAL_9 = SHELL_QUAD_9
//   , TET_4,  TETRAHEDRON_4  = TET_4
//   , TET_8,  TETRAHEDRON_8  = TET_8
//   , TET_10, TETRAHEDRON_10 = TET_10
//   , TET_11, TETRAHEDRON_11 = TET_11
//   , PYRAMID_5
//   , PYRAMID_13
//   , PYRAMID_14
//   , WEDGE_6
//   , WEDGE_12
//   , WEDGE_15
//   , WEDGE_18
//   , HEX_8,  HEXAHEDRON_8  = HEX_8
//   , HEX_20, HEXAHEDRON_20 = HEX_20
//   , HEX_27, HEXAHEDRON_27 = HEX_27
//   , END_TOPOLOGY
//   , NUM_TOPOLOGIES = END_TOPOLOGY - BEGIN_TOPOLOGY
//   , SUPEREDGE_START = END_TOPOLOGY+1
//   , SUPEREDGE_END = SUPEREDGE_START + 1000
//   , SUPERFACE_START = SUPEREDGE_END+1
//   , SUPERFACE_END = SUPERFACE_START + 1000
//   , SUPERELEMENT_START = SUPERFACE_END+1
//   , FORCE_TOPOLOGY_TO_UNSIGNED = ~0U // max unsigned int
// };

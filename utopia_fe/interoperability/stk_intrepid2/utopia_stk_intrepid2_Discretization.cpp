#include "utopia_stk_intrepid2_Discretization.hpp"

#include "utopia_stk_intrepid2.hpp"

namespace utopia {

    class Discretization<stk_FS_t, stk_FE_t>::Impl {
    public:
        std::shared_ptr<stk_FE_t> fe;

        // static ::shards::CellTopology convert_elem_type(::stk::topology::topology_t topo,
        //                                                 bool no_shell_topologies = false) {
        //     switch (topo) {
        //         case ::stk::topology::NODE:
        //             return ::shards::getCellTopologyData<::shards::Node>();

        //         case ::stk::topology::LINE_2_1D:
        //             return ::shards::getCellTopologyData<::shards::Line<>>();
        //         case ::stk::topology::BEAM_2:
        //         case ::stk::topology::LINE_2:
        //             return (no_shell_topologies ? ::shards::getCellTopologyData<::shards::Line<>>()
        //                                         : ::shards::getCellTopologyData<::shards::ShellLine<>>());
        //         case ::stk::topology::SHELL_LINE_2:
        //             return (no_shell_topologies ? ::shards::getCellTopologyData<::shards::Line<>>()
        //                                         : ::shards::getCellTopologyData<::shards::ShellLine<>>());
        //         case ::stk::topology::TRI_3_2D:
        //             return ::shards::getCellTopologyData<::shards::Triangle<>>();
        //         case ::stk::topology::TRI_3:
        //             return (no_shell_topologies ? ::shards::getCellTopologyData<::shards::Triangle<>>()
        //                                         : ::shards::getCellTopologyData<::shards::ShellTriangle<>>());
        //         case ::stk::topology::SHELL_TRI_3:
        //             return (no_shell_topologies ? ::shards::getCellTopologyData<::shards::Triangle<>>()
        //                                         : ::shards::getCellTopologyData<::shards::ShellTriangle<>>());
        //         case ::stk::topology::QUAD_4:
        //             return ::shards::getCellTopologyData<::shards::Quadrilateral<>>();
        //         case ::stk::topology::QUAD_4_2D:
        //             return ::shards::getCellTopologyData<::shards::Quadrilateral<>>();
        //         case ::stk::topology::SHELL_QUAD_4:
        //             return (no_shell_topologies ? ::shards::getCellTopologyData<::shards::Quadrilateral<>>()
        //                                         : ::shards::getCellTopologyData<::shards::ShellQuadrilateral<>>());
        //         case ::stk::topology::HEX_8:
        //             return ::shards::getCellTopologyData<::shards::Hexahedron<>>();
        //         case ::stk::topology::HEX_27:
        //             return ::shards::getCellTopologyData<::shards::Hexahedron<27>>();
        //         case ::stk::topology::TET_4:
        //             return ::shards::getCellTopologyData<::shards::Tetrahedron<>>();
        //         default: {
        //             assert(false);
        //             return ::shards::getCellTopologyData<::shards::Node>();
        //         }
        //     }
        // }

        // static bool create_from_buckets(const ::stk::mesh::BulkData &bulk_data,
        //                                 const BucketVector_t &buckets,
        //                                 FE &fe,
        //                                 int degree) {
        //     if (buckets.begin() == buckets.end()) {
        //         // utopia::err() << "[Warning] buckets.begin() == buckets.end()\n";
        //         return false;
        //     }

        //     auto &meta_data = bulk_data.mesh_meta_data();

        //     auto *first_bucket = *buckets.begin();

        //     // Dirty hack (FIXME once stk usage is a bit more profficient)
        //     auto topo = convert_elem_type(first_bucket->topology(), true);
        //     Size_t n_nodes_x_elem = bulk_data.num_nodes((*first_bucket)[0]);
        //     Size_t spatial_dim = meta_data.spatial_dimension();
        //     auto *coords = meta_data.coordinate_field();

        //     Size_t n_cells = 0;
        //     for (auto *b_ptr : buckets) {
        //         n_cells += b_ptr->size();
        //     }

        //     StkViewDevice_t<Scalar> device_cell_points("cell_points", n_cells, n_nodes_x_elem, spatial_dim);
        //     StkIntViewDevice_t device_element_tags("element_tags", n_cells);

        //     typename StkViewDevice_t<Scalar>::HostMirror cell_points =
        //     Kokkos::create_mirror_view(device_cell_points); StkIntViewDevice_t::HostMirror element_tags =
        //     Kokkos::create_mirror_view(device_element_tags);

        //     int elem_idx = 0;
        //     for (const auto &ib : buckets) {
        //         const Bucket_t &b = *ib;

        //         const Bucket_t::size_type length = b.size();

        //         const int block = utopia::stk::extract_set_id_from_bucket(b, b.topology().rank());

        //         for (Bucket_t::size_type k = 0; k < length; ++k) {
        //             Entity_t elem = b[k];
        //             const Size_t n_nodes = bulk_data.num_nodes(elem);

        //             element_tags(elem_idx) = block;

        //             assert(n_nodes == n_nodes_x_elem);

        //             auto node_ids = bulk_data.begin_nodes(elem);

        //             for (Size_t i = 0; i < n_nodes; ++i) {
        //                 const Scalar_t *point = (const Scalar_t *)::stk::mesh::field_data(*coords, node_ids[i]);

        //                 for (int d = 0; d < spatial_dim; ++d) {
        //                     cell_points(elem_idx, i, d) = point[d];
        //                 }
        //             }

        //             ++elem_idx;
        //         }
        //     }

        //     Kokkos::deep_copy(device_cell_points, cell_points);
        //     Kokkos::deep_copy(device_element_tags, element_tags);

        //     fe.init(topo, device_cell_points, degree);
        //     fe.element_tags() = device_element_tags;
        //     return true;
        // }
    };

    Discretization<stk_FS_t, stk_FE_t>::Discretization(const std::shared_ptr<FunctionSpace> &space,
                                                       const std::shared_ptr<FE> &fe)
        : Super(space), impl_(utopia::make_unique<Impl>()) {
        impl_->fe = fe;
    }

    Discretization<stk_FS_t, stk_FE_t>::~Discretization() {}

    void Discretization<stk_FS_t, stk_FE_t>::read(Input &in) {
        // TODO
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<stk_FS_t, stk_FE_t>::create(std::vector<FE_ptr> &fe, int order, const Part &part) {
        // FIXME
        FE_ptr fe0 = std::make_shared<FE>();
        create_fe(*this->space(), *fe0, order);

        fe.clear();
        fe.push_back(fe0);
    }

    void Discretization<stk_FS_t, stk_FE_t>::create_on_boundary(std::vector<FE_ptr> &fe, int order, const Part &part) {
        // FIXME
        FE_ptr fe0 = std::make_shared<FE>();

        if (part.is_all()) {
            create_fe_on_boundary(*this->space(), *fe0, order);
        } else {
            create_fe_on_boundary(*this->space(), *fe0, part.name, order);
        }

        fe.clear();
        fe.push_back(fe0);
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<stk_FS_t, stk_FE_t>::convert_field(const Field &in,
                                                           std::vector<std::shared_ptr<FEField>> &out,
                                                           const Part &part) {
        // FIXME

        if (out.size() == 1) {
            utopia::convert_field(in, *out[0]);
        } else {
            assert(false);
            // auto field0 = std::make_shared<FEField>();
            // utopia::convert_field(*field0, out);

            // out.clear();
            // out.push_back(field0);
        }
    }
    void Discretization<stk_FS_t, stk_FE_t>::convert_field(const std::vector<std::shared_ptr<FEField>> &in,
                                                           Field &out,
                                                           const Part &part) {
        // FIXME
        assert(in.size() == 1);
        utopia::convert_field(*in[0], out);
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<stk_FS_t, stk_FE_t>::global_to_local(const Vector &vector,
                                                             std::vector<VectorAccumulator> &element_vectors,
                                                             const Part &part,
                                                             const int comp) {
        // FIXME
        element_vectors.resize(1);
        utopia::global_to_local(*this->space(), vector, element_vectors[0], comp);
    }

    void Discretization<stk_FS_t, stk_FE_t>::local_to_global(const std::vector<MatrixAccumulator> &acc,
                                                             AssemblyMode mode,
                                                             Matrix &mat,
                                                             const Part &part) {
        // FIXME
        assert(acc.size() == 1);
        utopia::local_to_global(*this->space(), acc[0], mode, mat);
    }

    void Discretization<stk_FS_t, stk_FE_t>::local_to_global(const std::vector<VectorAccumulator> &acc,
                                                             AssemblyMode mode,
                                                             Vector &vec,
                                                             const Part &part) {
        // FIXME
        assert(acc.size() == 1);
        utopia::local_to_global(*this->space(), acc[0], mode, vec);
    }

    void Discretization<stk_FS_t, stk_FE_t>::local_to_global_on_boundary(const std::vector<VectorAccumulator> &acc,
                                                                         AssemblyMode mode,
                                                                         Vector &vec,
                                                                         const Part &part) {
        // TODO
        utopia::side_local_to_global(*this->space(), acc[0], mode, vec, part.name);
    }

}  // namespace utopia

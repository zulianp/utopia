#include "utopia_stk_DofMap.hpp"

#include "utopia_stk_Commons.hpp"

#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <cmath>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <vector>

std::string mpi_error_2_string(const int error_code) {
    std::string str;
    str.resize(MPI_MAX_ERROR_STRING);
    int len = MPI_MAX_ERROR_STRING;
    MPI_Error_string(error_code, &str[0], &len);
    str.resize(len);
    return str;
}

#define MPI_CATCH_ERROR(expr)                                                                                \
    {                                                                                                        \
        int ret = (expr);                                                                                    \
        if (MPI_SUCCESS != ret)                                                                              \
            std::cerr << "[" << utopia::mpi_world_rank() << "][Error] mpi error " << mpi_error_2_string(ret) \
                      << " code: " << ret << "(" << __FILE__ << ":" << __LINE__ << ")" << std::endl;         \
    }

namespace utopia {
    namespace stk {

        class DofMap::Impl {
        public:
            using SizeType = Traits<FunctionSpace>::SizeType;
            using IndexArray = Traits<FunctionSpace>::IndexArray;
            using Communicator = Traits<FunctionSpace>::Communicator;
            using Bucket = ::stk::mesh::Bucket;
            using BucketVector = ::stk::mesh::BucketVector;
            using Entity = ::stk::mesh::Entity;

            Communicator comm;
            IndexArray d_nnz, o_nnz;
            IndexArray local_to_global;
            SizeType owned_dof_start{-1}, owned_dof_end{-1};
            SizeType aura_nodes_offset{0};
            int n_var{1};
            bool valid_local_id_mode{false};

            void print_map(const ::stk::mesh::BulkData &bulk_data, std::ostream &os) const {
                auto &meta_data = bulk_data.mesh_meta_data();

                (!meta_data.locally_owned_part() & meta_data.globally_shared_part());

                // auto &node_buckets = utopia::stk::universal_nodes(bulk_data);

                os << "Owned\n";
                print_map(bulk_data, meta_data.locally_owned_part(), os);

                os << "Shared\n";
                print_map(bulk_data, (!meta_data.locally_owned_part() & meta_data.globally_shared_part()), os);

                os << "Aura\n";
                print_map(bulk_data, meta_data.aura_part(), os);
            }

            void print_map(const ::stk::mesh::BulkData &bulk_data,
                           const ::stk::mesh::Selector &selector,
                           std::ostream &os) const {
                const auto &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, selector);

                for (auto *b_ptr : node_buckets) {
                    auto &b = *b_ptr;
                    const SizeType length = b.size();

                    for (SizeType i = 0; i < length; ++i) {
                        auto node = b[i];

                        auto local_index = utopia::stk::convert_entity_to_index(node);

                        os << local_index << " -> " << bulk_data.local_id(node) << " ";

                        if (!local_to_global.empty()) {
                            auto dof = local_to_global[local_index];
                            os << dof << ", ";
                        }

                        os << utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node)) << '\n';
                    }
                }

                // os << "Linear local_to_global:\n";
                // if (!local_to_global.empty()) {
                //     SizeType n = local_to_global.size();
                //     for (SizeType i = 0; i < n; ++i) {
                //         os << i << " -> "
                //            << local_to_global[i]
                //            // << ", " << bulk_data.identifier(Entity(utopia::stk::convert_index_to_stk_index(i)))
                //            << "\n";
                //     }
                // }
            }
        };

        typedef struct GhostInfo {
            // public:
            using SizeType = Traits<FunctionSpace>::SizeType;
            SizeType g_id{-1};
            SizeType dof{-1};

            inline bool operator<(const GhostInfo &other) const { return g_id < other.g_id; }
        } GhostInfo;

        class MPIGhostInfo {
        public:
            MPIGhostInfo() {
                lengths[0] = 1;
                lengths[1] = 1;

                displacements[0] = static_cast<MPI_Aint>(offsetof(GhostInfo, g_id));
                displacements[1] = static_cast<MPI_Aint>(offsetof(GhostInfo, dof));

                types[0] = MPIType<GhostInfo::SizeType>::value();
                types[1] = MPIType<GhostInfo::SizeType>::value();

                MPI_CATCH_ERROR(MPI_Type_create_struct(2, lengths, displacements, types, &data_type));
                MPI_CATCH_ERROR(MPI_Type_commit(&data_type));

                // MPI_CATCH_ERROR(MPI_Type_contiguous(2, MPIType<GhostInfo::SizeType>::value(), &data_type));
                // MPI_CATCH_ERROR(MPI_Type_commit(&data_type));
            }

            ~MPIGhostInfo() { MPI_CATCH_ERROR(MPI_Type_free(&data_type)); }

            inline MPI_Datatype &get() { return data_type; }

        private:
            int lengths[2];
            MPI_Aint displacements[2];
            MPI_Datatype types[2];
            MPI_Datatype data_type;
        };

        class DofMap::DofExchange : public Describable {
        public:
            void increment_incoming(const int source) { ++rank_incoming[source]; }

            void increment_outgoing(const int destination) { ++rank_outgoing[destination]; }

            void add_dof_mapping_to_outgoing(const int destination,
                                             const Impl::SizeType &identifier,
                                             const Impl::SizeType &dof) {
                buffers_outgoing[destination].push_back({identifier, dof});
            }

            void add_global_to_o_nnz(const int destination,
                                     const Impl::SizeType &global_id,
                                     const Impl::SizeType &o_nnz) {
                buffers_outgoing[destination].push_back({global_id, o_nnz});
            }

            Impl::SizeType find_dof_from_incoming(const int source, const Impl::SizeType &identifier) {
                assert(source < int(buffers_incoming.size()));
                const auto &inc = buffers_incoming[source];
                assert(!inc.empty());

                if (inc.size() == 1) {
                    assert(inc[0].g_id == identifier);
                    return inc[0].dof;
                }

                GhostInfo wrapper;
                wrapper.g_id = identifier;

                const auto it = std::lower_bound(inc.begin(), inc.end(), wrapper);
#ifndef NDEBUG
                if (it->g_id != identifier) {
                    utopia::err() << "\n\n[" << comm.rank() << "] (source=" << source
                                  << ") Could not find dof for identifier: " << identifier << "!\n\n";
                }
#endif
                assert(it != inc.end());
                assert(it->g_id == identifier);

                return it->dof;
            }

            DofExchange(const Impl::Communicator &comm) : comm(comm) {}

            ~DofExchange() {}

            void describe(std::ostream &os) const {
                os << "outgoing:\n";
                for (auto p : rank_outgoing) {
                    os << p.first << " " << p.second << '\n';
                }

                os << "incoming:\n";
                for (auto p : rank_incoming) {
                    os << p.first << " " << p.second << '\n';
                }
            }

            void allocate() {
                buffers_outgoing.resize(comm.size());
                buffers_incoming.resize(comm.size());

                {
                    for (auto ro : rank_outgoing) {
                        assert(ro.second > 0);
                        assert(ro.first < int(buffers_outgoing.size()));

                        buffers_outgoing[ro.first].reserve(ro.second);
                    }
                }

                {
                    for (auto ri : rank_incoming) {
                        assert(ri.second > 0);
                        assert(ri.first < int(buffers_incoming.size()));
                        buffers_incoming[ri.first].resize(ri.second);
                    }
                }
            }

            void sort_outgoing() {
                for (auto &buff : buffers_outgoing) {
                    std::sort(std::begin(buff), std::end(buff));
                }
            }

            void exchange() {
                std::vector<MPI_Request> requests;
                requests.reserve(rank_incoming.size() + rank_outgoing.size());

                // auto raw_comm = MPI_COMM_WORLD;

                auto raw_comm = comm.raw_comm();

                for (auto ri : rank_incoming) {
                    requests.push_back(MPI_REQUEST_NULL);

                    assert(std::size_t(ri.first) < buffers_incoming.size());

                    auto &buff = buffers_incoming[ri.first];
                    assert(!buff.empty());

                    MPI_CATCH_ERROR(MPI_Irecv(&buff[0],
                                              buff.size(),
                                              data_type_.get(),
                                              ri.first,
                                              exchange_number_,
                                              raw_comm,
                                              &requests.back()));
                }

                for (auto ro : rank_outgoing) {
                    requests.push_back(MPI_REQUEST_NULL);

                    assert(std::size_t(ro.first) < buffers_outgoing.size());

                    auto &buff = buffers_outgoing[ro.first];
                    assert(!buff.empty());
                    assert(std::size_t(ro.second) == buff.size());

                    MPI_CATCH_ERROR(MPI_Isend(&buff[0],
                                              buff.size(),
                                              data_type_.get(),
                                              ro.first,
                                              exchange_number_,
                                              raw_comm,
                                              &requests.back()));
                }

                if (!requests.empty()) {
                    MPI_CATCH_ERROR(MPI_Waitall(static_cast<int>(requests.size()), &requests[0], MPI_STATUS_IGNORE));

                    ++exchange_number_;
                }
            }

            void describe_incoming(std::ostream &os) const {
                os << "in:\n";
                for (auto ri : rank_incoming) {
                    assert(ri.first < int(buffers_incoming.size()));
                    auto &buff = buffers_incoming[ri.first];

                    os << "rank: " << ri.first << " (size=" << buff.size() << "): ";
                    // for (auto b : buff) {
                    //     os << b.g_id << " -> " << b.dof << ", ";
                    // }

                    os << "\n";
                }
            }

            void describe_outgoing(std::ostream &os) const {
                os << "out:\n";
                for (auto ro : rank_outgoing) {
                    if (ro.first >= int(buffers_outgoing.size())) {
                        os << "Bad value: " << ro.first << " " << ro.second << '\n';
                        continue;
                    }

                    auto &buff = buffers_outgoing[ro.first];

                    os << "rank: " << ro.first << " (size=" << buff.size() << "): ";
                    // for (auto b : buff) {
                    //     os << b.g_id << " -> " << b.dof << ", ";
                    // }

                    os << "\n";
                }
            }

            void clear() {
                rank_incoming.clear();
                rank_outgoing.clear();
                clear_buffers();
            }

            void clear_counters() {
                for (auto &ri : rank_incoming) {
                    ri.second = 0;
                }

                for (auto &ro : rank_outgoing) {
                    ro.second = 0;
                }
            }

            void swap_incoming_outgoing() {
                std::swap(buffers_outgoing, buffers_incoming);
                std::swap(rank_incoming, rank_outgoing);

                for (auto &bo : buffers_outgoing) {
#ifndef NDEBUG
                    auto capacity = bo.capacity();
#endif  // NDEBUG

                    bo.resize(0);
                    assert(capacity == bo.capacity());
                }
            }

            void clear_buffers() {
                buffers_outgoing.clear();
                buffers_incoming.clear();
            }

            void add_to_o_nnz(const SizeType &local_offset, IndexArray &o_nnz) {
                for (auto &o : buffers_incoming) {
                    for (auto nnz : o) {
                        assert(nnz.g_id >= 0);
                        assert(nnz.g_id - local_offset >= 0);
                        assert(nnz.g_id - local_offset < SizeType(o_nnz.size()));

                        o_nnz[nnz.g_id - local_offset] += nnz.dof;
                    }
                }
            }

        private:
            const Impl::Communicator &comm;
            std::unordered_map<int, Impl::SizeType> rank_outgoing, rank_incoming;
            std::vector<std::vector<GhostInfo>> buffers_outgoing, buffers_incoming;

            MPIGhostInfo data_type_;
            int exchange_number_{0};
        };

        DofMap::DofMap() : impl_(utopia::make_unique<Impl>()) {}

        DofMap::~DofMap() = default;

        static void print_nodes(const ::stk::mesh::BulkData &bulk_data,
                                const ::stk::mesh::Selector &selector,
                                std::ostream &os) {
            using Bucket_t = ::stk::mesh::Bucket;
            using Entity_t = ::stk::mesh::Entity;

            const auto &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, selector);

            std::vector<int> procs;
            for (auto *b_ptr : node_buckets) {
                const auto &b = *b_ptr;
                const auto length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    const Entity_t &elem = b[k];
                    const auto local_id = utopia::stk::convert_entity_to_index(elem);
                    const auto id = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(elem));
                    os << "lid: " << local_id << ", id: " << id << ' ' << " => " << (bulk_data.local_id(elem))
                       << ", offset: " << elem.local_offset() << ", in range: " << bulk_data.in_index_range(elem)
                       << ", is_aura: " << bulk_data.is_aura_ghosted_onto_another_proc(bulk_data.entity_key(elem))
                       << ", is_shared: " << bulk_data.in_shared(elem) << ", owner: ";
                    os << bulk_data.parallel_owner_rank(elem);
                    os << ", recv: " << bulk_data.in_receive_ghost(elem) << ", send: " << bulk_data.in_send_ghost(elem);
                    os << '\n';

                    bulk_data.comm_procs(elem, procs);

                    os << "procs: ";
                    for (auto p : procs) {
                        os << p << ' ';
                    }

                    os << '\n';
                }

                os << '\n';
            }

            os << '\n';
        }

        void print_elements(const ::stk::mesh::BulkData &bulk_data,
                            const ::stk::mesh::Selector &selector,
                            std::ostream &os) {
            using Bucket_t = ::stk::mesh::Bucket;
            using Entity_t = ::stk::mesh::Entity;

            const auto &elem_buckets = bulk_data.get_buckets(::stk::topology::ELEMENT_RANK, selector);

            for (auto *b_ptr : elem_buckets) {
                const auto &b = *b_ptr;
                const auto length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    const Entity_t &elem = b[k];

                    const auto id = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(elem));
                    os << "id: " << id << ", owner: " << bulk_data.parallel_owner_rank(elem);
                    // os <<

                    os << '\n';

                    for (auto it = bulk_data.begin_nodes(elem); it < bulk_data.end_nodes(elem); ++it) {
                        os << '(' << utopia::stk::convert_stk_index_to_index(bulk_data.identifier(*it)) << ','
                           << bulk_data.parallel_owner_rank(*it) << ") ";
                    }

                    os << '\n';
                }
            }
        }

        int DofMap::n_var() const { return impl_->n_var; }
        void DofMap::set_n_var(const int n_var) { impl_->n_var = n_var; }

        void DofMap::set_valid_local_id_mode(const bool val) { impl_->valid_local_id_mode = val; }

        void DofMap::init(::stk::mesh::BulkData &bulk_data) {
            impl_->comm = Impl::Communicator(bulk_data.parallel());

            if (impl_->comm.size() == 1) {
                init_serial(bulk_data);
            } else {
                init_parallel(impl_->comm, bulk_data);
            }
        }

        void DofMap::init(Mesh &mesh, const bool verbose) {
            auto &meta_data = mesh.meta_data();
            auto &bulk_data = mesh.bulk_data();

            const int rank = mesh.comm().rank();
            const int size = mesh.comm().size();

            impl_->comm = Impl::Communicator(bulk_data.parallel());

            // std::stringstream ss;
            // mesh.describe(ss);

            // describe_mesh_connectivity(mesh);

            // TODO IF bulk_data.local_id(aura_node) != Crap then use it!

            if (mesh.has_aura() && size > 1) {
                // if (true) {
                UTOPIA_TRACE_REGION_BEGIN("DofMap::init_aura");

                using Bucket_t = ::stk::mesh::Bucket;
                using BucketVector_t = ::stk::mesh::BucketVector;
                using Entity_t = ::stk::mesh::Entity;

                auto &meta_data = bulk_data.mesh_meta_data();

                const ::stk::mesh::Selector universal_selector = meta_data.universal_part();
                const ::stk::mesh::Selector local_selector = meta_data.locally_owned_part();
                const BucketVector_t &elem_buckets =
                    bulk_data.get_buckets(::stk::topology::ELEMENT_RANK, local_selector);

                // print_nodes(bulk_data, universal_selector, ss);
                // print_elements(bulk_data, universal_selector, ss);

                SizeType n_local_nodes = count_local_nodes(bulk_data);
                SizeType n_universal_nodes = count_universal_nodes(bulk_data);

                SizeType offset = 0;
                mesh.comm().exscan_sum(&n_local_nodes, &offset, 1);

                std::vector<std::unordered_set<SizeType>> node2node(n_universal_nodes);

                for (auto *b_ptr : elem_buckets) {
                    const auto &b = *b_ptr;
                    const auto length = b.size();

                    for (Bucket_t::size_type k = 0; k < length; ++k) {
                        const Entity_t elem = b[k];
                        const auto node_ids = bulk_data.begin_nodes(elem);
                        const Size_t n_nodes = bulk_data.num_nodes(elem);

                        for (Size_t i = 0; i < n_nodes; ++i) {
                            const SizeType node_i = bulk_data.local_id(node_ids[i]);
                            assert(node_i < n_universal_nodes);
                            for (Size_t j = 0; j < n_nodes; ++j) {
                                node2node[node_i].insert(utopia::stk::convert_entity_to_index(node_ids[j]));
                            }
                        }
                    }
                }

                // mesh.comm().synched_print(ss.str());

                impl_->d_nnz.resize(n_local_nodes, 0);
                impl_->o_nnz.resize(n_local_nodes, 0);
                impl_->local_to_global.resize(n_universal_nodes, -1);
                impl_->owned_dof_start = offset;
                impl_->owned_dof_end = offset + n_local_nodes;
                //
                {
                    // Handle aura
                    const ::stk::mesh::Selector aura_selector = meta_data.aura_part();

                    const BucketVector_t &elem_buckets =
                        bulk_data.get_buckets(::stk::topology::ELEMENT_RANK, aura_selector);

                    for (auto *b_ptr : elem_buckets) {
                        const auto &b = *b_ptr;
                        const auto length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            const Entity_t elem = b[k];
                            const auto node_ids = bulk_data.begin_nodes(elem);
                            const Size_t n_nodes = bulk_data.num_nodes(elem);

                            for (Size_t i = 0; i < n_nodes; ++i) {
                                if (!bulk_data.in_receive_ghost(node_ids[i])) {
                                    const SizeType node_i = bulk_data.local_id(node_ids[i]);

                                    assert(node_i < n_universal_nodes);

                                    // if (node_i < n_universal_nodes) {
                                    for (Size_t j = 0; j < n_nodes; ++j) {
                                        node2node[node_i].insert(utopia::stk::convert_entity_to_index(node_ids[j]));

                                        // if (!bulk_data.in_receive_ghost(node_ids[j])) {
                                        // }
                                    }
                                }
                            }
                        }
                    }
                }

                {
                    // Non-zero pattern
                    const BucketVector_t &node_buckets =
                        bulk_data.get_buckets(::stk::topology::NODE_RANK, universal_selector);

                    Impl::SizeType max_local_node_id{0};
                    Impl::SizeType min_aura_node_id = std::numeric_limits<Impl::SizeType>::max();

                    Impl::SizeType index = 0;

                    for (auto *b_ptr : node_buckets) {
                        const auto &b = *b_ptr;
                        const auto length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            const Entity_t e = b[k];
                            if (bulk_data.in_receive_ghost(e)) {
                                min_aura_node_id =
                                    std::min(Impl::SizeType(utopia::stk::convert_entity_to_index(e)), min_aura_node_id);
                                continue;
                            }
                            Impl::SizeType i = bulk_data.local_id(e);
                            max_local_node_id = std::max(i, max_local_node_id);

                            int owner_rank = bulk_data.parallel_owner_rank(e);

                            auto &nodes = node2node[i];

                            if (owner_rank == rank) {
                                impl_->local_to_global[i] = offset + index;

                                for (auto n : nodes) {
                                    Impl::Entity node_e(utopia::stk::convert_index_to_stk_index(n));
                                    int node_owner_rank = bulk_data.parallel_owner_rank(node_e);

                                    if (node_owner_rank == owner_rank) {
                                        impl_->d_nnz[index] += 1;

                                    } else {
                                        impl_->o_nnz[index] += 1;
                                    }
                                }

                                ++index;
                            }
                        }
                    }

                    // ss << "max_local_node_id: " << max_local_node_id << "\n";
                    // ss << "min_aura_node_id: " << min_aura_node_id << "\n";

                    if (!impl_->valid_local_id_mode) {
                        impl_->aura_nodes_offset = (min_aura_node_id - (max_local_node_id + 1));
                        assert(impl_->aura_nodes_offset >= 0);
                    } else {
                        impl_->aura_nodes_offset = 0;
                    }

                    if (verbose) std::cout << "aura_nodes_offset: " << impl_->aura_nodes_offset << "\n";
                }

                exchange_shared_dofs(mesh.comm(), bulk_data, false);

                // mesh.describe(ss);
                // describe(ss);

                // mesh.comm().synched_print(ss.str());

                UTOPIA_TRACE_REGION_END("DofMap::init_aura");
            } else {
                init(bulk_data);
            }

            // mesh.describe(utopia::out().stream());
            // describe(ss);
            // describe_debug(mesh, ss);

            // mesh.comm().synched_print(ss.str());

            if (verbose) {
                std::stringstream ss;
                impl_->print_map(bulk_data, ss);
                mesh.comm().synched_print(ss.str(), utopia::out().stream());
            }
        }

        void DofMap::describe_mesh_connectivity(Mesh &mesh) {
            auto &meta_data = mesh.meta_data();
            auto &bulk_data = mesh.bulk_data();
            auto &comm = mesh.comm();
            const bool has_aura = bulk_data.is_automatic_aura_on();

            const int rank = comm.rank();
            const int size = comm.size();

            std::stringstream ss;

            auto &shared_node_buckets = shared_nodes(bulk_data);
            std::vector<SizeType> shared_out(size);
            std::vector<SizeType> shared_in(size);
            std::vector<SizeType> aura_out(size);
            std::vector<SizeType> aura_in(size);

            std::vector<int> procs;
            for (auto *b_ptr : shared_node_buckets) {
                const auto &b = *b_ptr;
                const auto length = b.size();

                for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                    const Impl::Entity node = b[k];
                    int owner_rank = bulk_data.parallel_owner_rank(node);

                    bulk_data.comm_procs(node, procs);

                    if (owner_rank == rank) {
                        for (auto p : procs) {
                            if (!bulk_data.in_send_ghost(bulk_data.entity_key(node), p)) {
                                ++shared_out[p];
                            }
                        }
                    } else {
                        ++shared_in[owner_rank];
                    }
                }
            }

            ss << "shared_out\n";

            for (int i = 0; i < size; ++i) {
                ss << shared_out[i];

                if (i < size - 1) {
                    ss << ',';
                }
            }

            ss << '\n';

            ss << "shared_in\n";

            for (int i = 0; i < size; ++i) {
                ss << shared_in[i];

                if (i < size - 1) {
                    ss << ',';
                }
            }

            ss << '\n';

            if (has_aura) {
                for (auto *b_ptr : universal_nodes(bulk_data)) {
                    const auto &b = *b_ptr;
                    const auto length = b.size();

                    for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                        const Impl::Entity node = b[k];

                        int owner_rank = bulk_data.parallel_owner_rank(node);

                        bulk_data.comm_procs(node, procs);

                        if (bulk_data.in_send_ghost(node)) {
                            for (auto p : procs) {
                                if (bulk_data.in_send_ghost(bulk_data.entity_key(node), p)) ++aura_out[p];
                            }

                        } else if (bulk_data.in_receive_ghost(node)) {
                            ++aura_in[owner_rank];
                        }
                    }
                }
            }

            ss << "aura_out\n";

            for (int i = 0; i < size; ++i) {
                ss << aura_out[i];

                if (i < size - 1) {
                    ss << ',';
                }
            }

            ss << '\n';

            ss << "aura_in\n";

            for (int i = 0; i < size; ++i) {
                ss << aura_in[i];

                if (i < size - 1) {
                    ss << ',';
                }
            }

            ss << '\n';

            comm.synched_print(ss.str());
        }

        void DofMap::exchange_shared_dofs(const Communicator &comm,
                                          ::stk::mesh::BulkData &bulk_data,
                                          bool add_to_nnz_pattern) {
            const int rank = comm.rank();
            const bool has_aura = bulk_data.is_automatic_aura_on();

            // auto &meta_data = bulk_data.mesh_meta_data();
            SizeType n_local_nodes = count_local_nodes(bulk_data);
            SizeType n_universal_nodes = count_universal_nodes(bulk_data);

            SizeType offset = 0;
            comm.exscan_sum(&n_local_nodes, &offset, 1);
            std::vector<std::unordered_set<SizeType>> node2node(n_universal_nodes);

            DofExchange dof_exchange(comm);

            auto &shared_node_buckets = shared_nodes(bulk_data);

            std::vector<int> procs;

            if (has_aura) {
                for (auto *b_ptr : universal_nodes(bulk_data)) {
                    const auto &b = *b_ptr;
                    const auto length = b.size();

                    for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                        const Impl::Entity node = b[k];

                        int owner_rank = bulk_data.parallel_owner_rank(node);

                        bulk_data.comm_procs(node, procs);

                        if (bulk_data.in_send_ghost(node)) {
                            for (auto p : procs) {
                                if (bulk_data.in_send_ghost(bulk_data.entity_key(node), p)) {
                                    dof_exchange.increment_outgoing(p);
                                }
                            }
                        } else if (bulk_data.in_receive_ghost(node)) {
                            dof_exchange.increment_incoming(owner_rank);
                        }
                    }
                }

                for (auto *b_ptr : shared_node_buckets) {
                    const auto &b = *b_ptr;
                    const auto length = b.size();

                    for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                        const Impl::Entity node = b[k];
                        int owner_rank = bulk_data.parallel_owner_rank(node);

                        bulk_data.comm_procs(node, procs);

                        if (owner_rank == rank) {
                            for (auto p : procs) {
                                if (!bulk_data.in_send_ghost(bulk_data.entity_key(node), p)) {
                                    dof_exchange.increment_outgoing(p);
                                }
                            }
                        } else {
                            dof_exchange.increment_incoming(owner_rank);
                        }
                    }
                }

            } else {
                for (auto *b_ptr : shared_node_buckets) {
                    const auto &b = *b_ptr;
                    const auto length = b.size();

                    for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                        const Impl::Entity node = b[k];
                        int owner_rank = bulk_data.parallel_owner_rank(node);

                        bulk_data.comm_procs(node, procs);

                        if (owner_rank == rank) {
                            for (auto p : procs) {
                                dof_exchange.increment_outgoing(p);
                            }
                        } else {
                            dof_exchange.increment_incoming(owner_rank);
                        }
                    }
                }
            }

            dof_exchange.allocate();

            if (has_aura) {
                for (auto *b_ptr : universal_nodes(bulk_data)) {
                    const auto &b = *b_ptr;
                    const auto length = b.size();

                    for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                        const Impl::Entity node = b[k];
                        int owner_rank = bulk_data.parallel_owner_rank(node);
                        auto local_index = utopia::stk::convert_entity_to_index(node);

                        if (bulk_data.in_send_ghost(node)) {
                            auto dof = impl_->local_to_global[local_index];
                            assert(dof >= 0);
                            bulk_data.comm_procs(node, procs);
                            for (auto p : procs) {
                                if (bulk_data.in_send_ghost(bulk_data.entity_key(node), p)) {
                                    dof_exchange.add_dof_mapping_to_outgoing(
                                        p, utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node)), dof);
                                }
                            }
                        }
                    }
                }

                for (auto *b_ptr : shared_node_buckets) {
                    const auto &b = *b_ptr;
                    const auto length = b.size();

                    for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                        const Impl::Entity node = b[k];
                        int owner_rank = bulk_data.parallel_owner_rank(node);

                        bulk_data.comm_procs(node, procs);

                        if (owner_rank == rank) {
                            auto local_index = utopia::stk::convert_entity_to_index(node);
                            auto dof = impl_->local_to_global[local_index];
                            assert(dof >= 0);

                            for (auto p : procs) {
                                if (!bulk_data.in_send_ghost(bulk_data.entity_key(node), p)) {
                                    dof_exchange.add_dof_mapping_to_outgoing(
                                        p, utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node)), dof);
                                }
                            }
                        }
                    }
                }

            } else {
                for (auto *b_ptr : shared_node_buckets) {
                    const auto &b = *b_ptr;
                    const auto length = b.size();

                    for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                        const Impl::Entity node = b[k];
                        int owner_rank = bulk_data.parallel_owner_rank(node);
                        auto local_index = utopia::stk::convert_entity_to_index(node);

                        if (owner_rank == rank) {
                            auto dof = impl_->local_to_global[local_index];
                            assert(dof >= 0);
                            bulk_data.comm_procs(node, procs);
                            for (auto p : procs) {
                                dof_exchange.add_dof_mapping_to_outgoing(
                                    p, utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node)), dof);
                            }
                        }
                    }
                }
            }

            dof_exchange.sort_outgoing();

            // std::stringstream ss;
            // dof_exchange.describe_outgoing(ss);

            // dof_exchange.describe_incoming(ss);

            // comm.synched_print(ss.str());

            dof_exchange.exchange();

            // dof_exchange.describe_incoming(ss);

            // comm.synched_print(ss.str());

            // comm.barrier();

            if (has_aura) {
                for (auto *b_ptr : aura_nodes(bulk_data)) {
                    const auto &b = *b_ptr;
                    const auto length = b.size();

                    for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                        const Impl::Entity node = b[k];
                        int owner_rank = bulk_data.parallel_owner_rank(node);

                        if (owner_rank != rank) {
                            auto g_id = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node));
                            auto dof = dof_exchange.find_dof_from_incoming(owner_rank, g_id);

                            auto local_index = utopia::stk::convert_entity_to_index(node) - impl_->aura_nodes_offset;

                            // assert(g_id == dof);
                            assert(dof >= 0);
                            assert(local_index >= 0);

                            if (!(local_index < impl_->local_to_global.size())) {
                                utopia::err() << local_index << "<" << impl_->local_to_global.size() << "\n";
                                assert(local_index < impl_->local_to_global.size());
                            }

                            impl_->local_to_global[local_index] = dof;
                        }
                    }
                }
            }
            // else {
            for (auto *b_ptr : shared_node_buckets) {
                const auto &b = *b_ptr;
                const auto length = b.size();

                for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                    const Impl::Entity node = b[k];
                    int owner_rank = bulk_data.parallel_owner_rank(node);

                    if (owner_rank != rank) {
                        auto g_id = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node));
                        auto dof = dof_exchange.find_dof_from_incoming(owner_rank, g_id);

                        auto local_index = utopia::stk::convert_entity_to_index(node);

                        // assert(g_id == dof);
                        assert(dof >= 0);
                        assert(local_index >= 0);
                        assert(local_index < impl_->local_to_global.size());
                        impl_->local_to_global[local_index] = dof;
                    }
                }
            }
            // }

            if (add_to_nnz_pattern) {
                dof_exchange.swap_incoming_outgoing();

                for (Impl::SizeType i = 0; i < n_universal_nodes; ++i) {
                    Impl::Entity e(utopia::stk::convert_index_to_stk_index(i));
                    if (bulk_data.in_receive_ghost(e)) continue;

                    int owner_rank = bulk_data.parallel_owner_rank(e);

                    if (owner_rank != rank) {
                        Impl::SizeType count = 0;
                        auto &nodes = node2node[i];

                        for (auto n : nodes) {
                            Impl::Entity node_e(utopia::stk::convert_index_to_stk_index(n));
                            int node_owner_rank = bulk_data.parallel_owner_rank(node_e);

                            if (node_owner_rank == rank) {
                                ++count;
                            }
                        }

                        assert(impl_->local_to_global[i] >= 0);
                        assert(impl_->local_to_global[i] < impl_->owned_dof_end);

                        dof_exchange.add_global_to_o_nnz(owner_rank, impl_->local_to_global[i], count);
                    }
                }

                comm.barrier();
                // {
                //     std::stringstream ss;
                //     // dof_exchange.describe(ss);
                //     dof_exchange.describe_outgoing(ss);
                //     comm.synched_print(ss.str());
                // }

                dof_exchange.exchange();
                dof_exchange.add_to_o_nnz(offset, impl_->o_nnz);
            }
        }  // namespace stk

        void DofMap::init_parallel(const Communicator &comm, ::stk::mesh::BulkData &bulk_data) {
            UTOPIA_TRACE_REGION_BEGIN("DofMap::init_parallel");
            const int rank = comm.rank();

            // auto &meta_data = bulk_data.mesh_meta_data();
            SizeType n_local_nodes = count_local_nodes(bulk_data);
            SizeType n_universal_nodes = count_universal_nodes(bulk_data);

            SizeType offset = 0;
            comm.exscan_sum(&n_local_nodes, &offset, 1);
            std::vector<std::unordered_set<SizeType>> node2node(n_universal_nodes);

            // auto &element_buckets = universal_elements(bulk_data);
            auto &element_buckets = local_elements(bulk_data);

            for (auto *b_ptr : element_buckets) {
                const auto &b = *b_ptr;
                const auto length = b.size();

                for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                    const Impl::Entity elem = b[k];
                    const auto node_ids = bulk_data.begin_nodes(elem);
                    const Size_t n_nodes = bulk_data.num_nodes(elem);

                    for (Size_t i = 0; i < n_nodes; ++i) {
                        const Impl::SizeType node_i =
                            utopia::stk::convert_stk_index_to_index(node_ids[i].local_offset());

                        for (Size_t j = 0; j < n_nodes; ++j) {
                            node2node[node_i].insert(
                                utopia::stk::convert_stk_index_to_index(node_ids[j].local_offset()));
                        }
                    }
                }
            }

            impl_->d_nnz.resize(n_local_nodes, 0);
            impl_->o_nnz.resize(n_local_nodes, 0);
            impl_->local_to_global.resize(n_universal_nodes, -1);
            impl_->owned_dof_start = offset;
            impl_->owned_dof_end = offset + n_local_nodes;

            Impl::SizeType index = 0;
            for (Impl::SizeType i = 0; i < n_universal_nodes; ++i) {
                auto &nodes = node2node[i];
                Impl::Entity e(utopia::stk::convert_index_to_stk_index(i));
                int owner_rank = bulk_data.parallel_owner_rank(e);

                if (owner_rank == rank) {
                    impl_->local_to_global[i] = offset + index;

                    for (auto n : nodes) {
                        Impl::Entity node_e(utopia::stk::convert_index_to_stk_index(n));
                        int node_owner_rank = bulk_data.parallel_owner_rank(node_e);

                        if (node_owner_rank == owner_rank) {
                            impl_->d_nnz[index] += 1;

                        } else {
                            impl_->o_nnz[index] += 1;
                        }
                    }

                    ++index;
                }
            }

            exchange_shared_dofs(comm, bulk_data);

            UTOPIA_TRACE_REGION_END("DofMap::init_parallel");
        }

        void DofMap::init_serial(::stk::mesh::BulkData &bulk_data) {
            UTOPIA_TRACE_REGION_BEGIN("DofMap::init_serial");
            using Bucket_t = ::stk::mesh::Bucket;
            using BucketVector_t = ::stk::mesh::BucketVector;
            using Entity_t = ::stk::mesh::Entity;

            auto &meta_data = bulk_data.mesh_meta_data();

            const ::stk::mesh::Selector selector = meta_data.locally_owned_part();
            const BucketVector_t &node_buckets = bulk_data.get_buckets(::stk::topology::ELEMENT_RANK, selector);

            const SizeType nln = count_local_nodes(bulk_data);
            std::vector<std::unordered_set<SizeType>> node2node(nln);

            for (auto *b_ptr : node_buckets) {
                const auto &b = *b_ptr;
                const auto length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    const Entity_t elem = b[k];
                    const auto node_ids = bulk_data.begin_nodes(elem);
                    const Size_t n_nodes = bulk_data.num_nodes(elem);

                    for (Size_t i = 0; i < n_nodes; ++i) {
                        // const auto node_i =
                        // utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node_ids[i]));
                        const auto node_i = utopia::stk::convert_entity_to_index(node_ids[i]);

                        for (Size_t j = 0; j < n_nodes; ++j) {
                            node2node[node_i].insert(
                                // utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node_ids[j])));
                                utopia::stk::convert_entity_to_index(node_ids[j]));
                        }
                    }
                }
            }

            impl_->d_nnz.resize(nln, 0);
            impl_->o_nnz.resize(nln, 0);

            for (SizeType i = 0; i < nln; ++i) {
                impl_->d_nnz[i] = node2node[i].size();
            }

            UTOPIA_TRACE_REGION_END("DofMap::init_serial");
        }

        bool DofMap::empty() const { return impl_->local_to_global.empty() && impl_->d_nnz.empty(); }

        void DofMap::describe(std::ostream &os) const {
            auto &d_nnz = impl_->d_nnz;
            auto &o_nnz = impl_->o_nnz;
            auto &local_to_global = impl_->local_to_global;

            Impl::SizeType index = 0;
            for (auto ltog : local_to_global) {
                os << index << " -> " << ltog << '\n';
                index++;
            }

            index = 0;
            for (auto d : d_nnz) {
                os << (impl_->owned_dof_start + index) << " -> nnz: " << d << ',' << o_nnz[index] << '\n';
                index++;
            }
        }

        void DofMap::describe_debug(Mesh &mesh, std::ostream &os) const {
            auto &d_nnz = impl_->d_nnz;
            auto &o_nnz = impl_->o_nnz;
            auto &local_to_global = impl_->local_to_global;

            Impl::SizeType index = 0;
            for (auto ltog : local_to_global) {
                if (ltog == -1) {
                    Impl::Entity e(utopia::stk::convert_index_to_stk_index(index));
                    // os << mesh.bulk_data().in_receive_ghost(e) << "\n";

                    os << index << " -> " << mesh.bulk_data().identifier(e) << '\t';
                    os << (mesh.bulk_data().in_shared(e) ? "shared " : "")
                       << "| owner: " << mesh.bulk_data().parallel_owner_rank(e) << "\n";
                }
                index++;
            }

            // index = 0;
            // for (auto d : d_nnz) {
            //     os << (impl_->owned_dof_start + index) << " -> nnz: " << d << ',' << o_nnz[index] << '\n';
            //     index++;
            // }
        }

        DofMap::GlobalIndex DofMap::local_to_global() const { return GlobalIndex(impl_->local_to_global, n_var()); }
        const DofMap::IndexArray &DofMap::d_nnz() const { return impl_->d_nnz; }
        const DofMap::IndexArray &DofMap::o_nnz() const { return impl_->o_nnz; }

        void DofMap::global_to_local(const Vector &global, Vector &local) const {
            if (impl_->comm.size() == 1) {
                local = global;
                return;
            }

            auto &l2g = impl_->local_to_global;

            ////////////////////////////////////////////////////////////
            SizeType n_owned_nodes = impl_->owned_dof_end - impl_->owned_dof_start;
            SizeType n_var_vec = global.local_size() / n_owned_nodes;

            assert(n_owned_nodes * n_var_vec == global.local_size());

            if (n_var_vec == 1) {
                local.zeros(layout(Communicator::self(), l2g.size(), l2g.size()));
                global.select(l2g, local);
            } else {
                local.zeros(layout(Communicator::self(), l2g.size() * n_var_vec, l2g.size() * n_var_vec));
                global.blocked_select(l2g, local, n_var_vec);
            }
            ////////////////////////////////////////////////////////////

            // if (n_var() == 1) {
            //     local.zeros(layout(Communicator::self(), l2g.size(), l2g.size()));
            //     global.select(l2g, local);
            // } else {
            //     const int nv = this->n_var();
            //     local.zeros(layout(Communicator::self(), l2g.size() * nv, l2g.size() * nv));
            //     global.blocked_select(l2g, local, nv);
            // }
        }

        void DofMap::local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const {
            assert(!global.empty());

            if (global.comm().size() == 1) {
                // assert(mpi_world_size() == 1);
                global = local;
                return;
            }

            switch (mode) {
                case OVERWRITE_MODE: {
                    auto l2g = this->local_to_global();

                    auto l_view = const_local_view_device(local);
                    auto g_view = local_view_device(global);

                    // TODO
                    // parallel_for(local_range_device(local), UTOPIA_LAMBDA(const SizeType){});

                    auto r = range(global);
                    const SizeType nl = l2g.size();
                    // const int nv = n_var();
                    const int nv = local.local_size() / nl;

                    assert(nv > 0);
                    assert(nl > 0);

                    assert(nv * nl == local.local_size());

                    for (SizeType i = 0; i < nl; ++i) {
                        for (int d = 0; d < nv; ++d) {
                            const SizeType g_id = l2g(i, d);

                            if (r.inside(g_id)) {
                                const SizeType l_id = i * nv + d;
                                assert(l_id < nv * nl);
                                assert(l_id < local.local_size());

                                g_view.set(g_id - r.begin(), l_view.get(l_id));
                            }
                        }
                    }

                    break;
                }
                default: {
                    assert(false && "IMPLEMENT ME");
                    Utopia::Abort();
                }
            }
        }

        DofMap::SizeType DofMap::shift_aura_idx(const SizeType idx) const { return idx - impl_->aura_nodes_offset; }

    }  // namespace stk
}  // namespace utopia

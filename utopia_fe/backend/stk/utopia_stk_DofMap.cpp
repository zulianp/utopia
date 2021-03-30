#include "utopia_stk_DofMap.hpp"

#include "utopia_stk_Commons.hpp"

#include "utopia_stk_FunctionSpace.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

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

            IndexArray d_nnz, o_nnz;
            IndexArray local_to_global;
            SizeType owned_dof_start{-1}, owned_dof_end{-1};
        };

        typedef struct IdToDof {
            // public:
            using SizeType = Traits<FunctionSpace>::SizeType;
            SizeType g_id{-1};
            SizeType dof{-1};

            inline bool operator<(const IdToDof &other) const { return g_id < other.g_id; }
        } IdToDof;

        class MPIIdToDof {
        public:
            MPIIdToDof() {
                lengths[0] = 1;
                lengths[1] = 1;

                displacements[0] = static_cast<MPI_Aint>(offsetof(IdToDof, g_id));
                displacements[1] = static_cast<MPI_Aint>(offsetof(IdToDof, dof));

                types[0] = MPIType<IdToDof::SizeType>::value();
                types[1] = MPIType<IdToDof::SizeType>::value();

                MPI_CATCH_ERROR(MPI_Type_create_struct(2, lengths, displacements, types, &data_type));
                MPI_CATCH_ERROR(MPI_Type_commit(&data_type));

                // MPI_CATCH_ERROR(MPI_Type_contiguous(2, MPIType<IdToDof::SizeType>::value(), &data_type));
                // MPI_CATCH_ERROR(MPI_Type_commit(&data_type));
            }

            ~MPIIdToDof() { MPI_CATCH_ERROR(MPI_Type_free(&data_type)); }

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

            Impl::SizeType find_dof_from_incoming(const int source, const Impl::SizeType &identifier) {
                assert(source < int(buffers_incoming.size()));
                const auto &inc = buffers_incoming[source];
                assert(!inc.empty());

                if (inc.size() == 1) {
                    assert(inc[0].g_id == identifier);
                    return inc[0].dof;
                }

                IdToDof wrapper;
                wrapper.g_id = identifier;

                const auto it = std::lower_bound(inc.begin(), inc.end(), wrapper);
                assert(it != inc.end());

                return it->dof;
            }

            DofExchange(const Impl::Communicator &comm) : comm(comm) {}

            ~DofExchange() {}

            void describe(std::ostream &os) const {
                os << "ougoing:\n";
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

                    auto &buff = buffers_incoming[ri.first];
                    assert(!buff.empty());

                    MPI_CATCH_ERROR(MPI_Irecv(
                        &buff[0], buff.size(), data_type_.get(), ri.first, comm.rank(), raw_comm, &requests.back()));
                }

                for (auto ro : rank_outgoing) {
                    requests.push_back(MPI_REQUEST_NULL);

                    auto &buff = buffers_outgoing[ro.first];
                    assert(!buff.empty());

                    MPI_CATCH_ERROR(MPI_Isend(
                        &buff[0], buff.size(), data_type_.get(), ro.first, ro.first, raw_comm, &requests.back()));
                }

                MPI_CATCH_ERROR(MPI_Waitall(static_cast<int>(requests.size()), &requests[0], MPI_STATUS_IGNORE));
            }

            void describe_incoming(std::ostream &os) const {
                for (auto ri : rank_incoming) {
                    assert(ri.first < int(buffers_incoming.size()));
                    auto &buff = buffers_incoming[ri.first];

                    os << ri.first << ": ";
                    for (auto b : buff) {
                        os << b.g_id << " -> " << b.dof << ", ";
                    }

                    os << "\n";
                }
            }

            void clear() {
                rank_outgoing.clear();
                rank_incoming.clear();
                buffers_outgoing.clear();
                buffers_incoming.clear();
            }

        private:
            const Impl::Communicator &comm;
            std::unordered_map<int, Impl::SizeType> rank_outgoing, rank_incoming;
            std::vector<std::vector<IdToDof>> buffers_outgoing, buffers_incoming;

            MPIIdToDof data_type_;
        };

        DofMap::DofMap() : impl_(utopia::make_unique<Impl>()) {}

        DofMap::~DofMap() = default;

        void print_nodes(const ::stk::mesh::BulkData &bulk_data,
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
                    const auto id = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(elem));
                    os << id << ' '
                       << bulk_data.parallel_owner_rank(elem)
                       // << ' ' << bulk_data.count_valid_connectivity(elem, ::stk::topology::NODE_RANK)
                       << '\n';

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
                    os << id << ' ' << bulk_data.parallel_owner_rank(elem) << '\n';

                    for (auto it = bulk_data.begin_nodes(elem); it < bulk_data.end_nodes(elem); ++it) {
                        os << '(' << utopia::stk::convert_stk_index_to_index(bulk_data.identifier(*it)) << ','
                           << bulk_data.parallel_owner_rank(*it) << ") ";
                    }

                    os << '\n';
                }
            }
        }

        void DofMap::init(::stk::mesh::BulkData &bulk_data) {
            Impl::Communicator comm(bulk_data.parallel());
            int rank = comm.rank();

            if (comm.size() == 1) return;

            auto &meta_data = bulk_data.mesh_meta_data();
            SizeType n_local_nodes = count_local_nodes(bulk_data);
            SizeType n_universal_nodes = count_universal_nodes(bulk_data);

            SizeType offset = 0;
            comm.exscan_sum(&n_local_nodes, &offset, 1);

            {
                std::stringstream ss;
                print_nodes(bulk_data, meta_data.globally_shared_part(), ss);
                comm.synched_print(ss.str());
            }

            std::vector<std::unordered_set<SizeType>> node2node(n_universal_nodes);

            auto &element_buckets = universal_elements(bulk_data);

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
                    for (auto n : nodes) {
                        Impl::Entity node_e(utopia::stk::convert_index_to_stk_index(n));
                        int node_owner_rank = bulk_data.parallel_owner_rank(node_e);

                        if (node_owner_rank == owner_rank) {
                            impl_->d_nnz[index] += 1;
                            impl_->local_to_global[index] = offset + index;
                        } else {
                            impl_->o_nnz[index] += 1;
                        }
                    }

                    ++index;
                } else {
                    // add o_nnz
                }
            }

            DofExchange dof_exchange(comm);

            auto &shared_node_buckets = shared_nodes(bulk_data);

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
                            dof_exchange.increment_outgoing(p);
                        }
                    } else {
                        for (auto p : procs) {
                            dof_exchange.increment_incoming(p);
                        }
                    }
                }
            }

            dof_exchange.allocate();

            index = 0;
            for (auto *b_ptr : shared_node_buckets) {
                const auto &b = *b_ptr;
                const auto length = b.size();

                for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                    const Impl::Entity node = b[k];
                    int owner_rank = bulk_data.parallel_owner_rank(node);

                    bulk_data.comm_procs(node, procs);

                    if (owner_rank == rank) {
                        for (auto p : procs) {
                            dof_exchange.add_dof_mapping_to_outgoing(
                                p, utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node)), offset + index);
                        }

                        index++;
                    }
                }
            }

            dof_exchange.sort_outgoing();
            dof_exchange.exchange();

            for (auto *b_ptr : shared_node_buckets) {
                const auto &b = *b_ptr;
                const auto length = b.size();

                for (Impl::Bucket::size_type k = 0; k < length; ++k) {
                    const Impl::Entity node = b[k];
                    int owner_rank = bulk_data.parallel_owner_rank(node);

                    if (owner_rank != rank) {
                        auto dof = dof_exchange.find_dof_from_incoming(
                            owner_rank, utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node)));

                        auto index = utopia::stk::convert_entity_to_index(node);
                        impl_->local_to_global[index] = dof;
                    }
                }
            }

            // TODO Comunicate o_nnz

            std::stringstream ss;
            // describe(ss);
            // dof_exchange.describe_incoming(ss);

            for (auto &nodes : node2node) {
                for (auto &n : nodes) {
                    ss << n << " ";
                }

                ss << "\n";
            }

            describe(ss);
            comm.synched_print(ss.str());

            Utopia::Abort();
        }

        bool DofMap::empty() const { return impl_->local_to_global.empty(); }

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
                os << (impl_->owned_dof_start + index) << " -> nnz: " << d << ", " << o_nnz[index] << '\n';
                index++;
            }
        }
    }  // namespace stk
}  // namespace utopia

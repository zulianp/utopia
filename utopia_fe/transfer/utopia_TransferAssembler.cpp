#include "utopia_TransferAssembler.hpp"

#include "utopia_LocalAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_QMortarBuilder.hpp"

#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/elem.h"

#include "utopia_libmesh.hpp"
#include "utopia_VTree.hpp"

#include "utopia_ElementDofMap.hpp"
#include "utopia_FESpacesAdapter.hpp"

#include "moonolith_profiler.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_tree.hpp"
#include "moonolith_n_tree_mutator_factory.hpp"
#include "moonolith_n_tree_with_span_mutator_factory.hpp"
#include "moonolith_n_tree_with_tags_mutator_factory.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "par_moonolith.hpp"

#include "utopia_Socket.hpp"

#include "utopia_private_PetrovGalerkinAssembler.hpp"

#include <cmath>
#include <queue>
#include <algorithm>
#include <sstream>
#include <numeric>



/**
 * TODO:
   - Allow element-node-dof <-> node-dof convertions
   - Construction of accumulation operators
   - Discrimination of element-matrices and deletion of bad entries (e.g., duplicate volume-surface maps)
 */

namespace utopia {

    template<int Dimensions>
    class FESpaceSerializerDeserializer {
    public:
        using InputStream  = moonolith::InputStream;
        using OutputStream = moonolith::OutputStream;

        using NTreeT 		= utopia::VTree<Dimensions>;
        using DataContainer = typename NTreeT::DataContainer;
        using Adapter       = typename NTreeT::DataType;

        FESpaceSerializerDeserializer(
            const libMesh::Parallel::Communicator &comm,
            const TransferOptions &opts,
            const std::shared_ptr<FESpacesAdapter> &local_spaces)
        : comm(comm),
          m_comm(comm.get()),
          opts(opts),
          local_spaces(local_spaces)
        {}

        const libMesh::Parallel::Communicator &comm;
        moonolith::Communicator m_comm;
        const TransferOptions &opts;

        std::shared_ptr<FESpacesAdapter> local_spaces;
        std::map<long, std::shared_ptr<FESpacesAdapter> > spaces;
        std::map<long, std::vector<std::shared_ptr<FESpacesAdapter> > > migrated_spaces;

        void read(
            const long ownerrank,
            const long senderrank,
            bool is_forwarding, DataContainer &data,
            InputStream &in
        ) {

            CHECK_STREAM_READ_BEGIN("vol_proj", in);

            std::shared_ptr<FESpacesAdapter> proc_space = std::make_shared<FESpacesAdapter>(m_comm);

            read_spaces(in, *proc_space, comm, comm);

            if (!is_forwarding) {
                assert(!spaces[ownerrank]);
                spaces[ownerrank] = proc_space;
            } else {
                migrated_spaces[ownerrank].push_back(proc_space);
            }

            data.reserve(data.size() + 3000);

            long offset = 0;

            if(opts.tags.empty()){
                int space_num = 0;
                for(auto s : proc_space->spaces()) {
                    if(s) {
                        //ID_FIX this should be fine n_elem is actually local sence the mesh is a SerialMesh
                        for (int i=0; i<s->n_elem(); i++) {
                            data.push_back(Adapter(*s, i, offset + i,space_num));
                            assert(!proc_space->dof_map(space_num)[i].empty());
                            data.back().set_dof_map(&proc_space->dof_map(space_num)[i].global);
                        }

                        offset += s->n_elem();

                    }

                    ++space_num;
                }
            } else {
                int space_num = 0;
                for(auto s : proc_space->spaces()) {
                    if(s) {
                        for (int i=0; i<s->n_elem(); i++) {
                            const libMesh::Elem * elem = s->elem_ptr(i);
                            //Volume Tag
                            int volume_tag = elem->subdomain_id();
                            data.push_back(Adapter(*s, i, offset + i,volume_tag) );
                            assert(!proc_space->dof_map(space_num)[i].empty());
                            data.back().set_dof_map(&proc_space->dof_map(space_num)[i].global);
                        }

                        offset += s->n_elem();

                    }

                    ++space_num;
                }
            }

            CHECK_STREAM_READ_END("vol_proj", in);
        };

        void write(
            const long ownerrank, const long recvrank,
            const std::vector<long>::const_iterator &begin,
            const std::vector<long>::const_iterator &end,
            const DataContainer &data,
            OutputStream &out) {

            CHECK_STREAM_WRITE_BEGIN("vol_proj", out);

            if (ownerrank == m_comm.rank()) {
                write_element_selection(begin, end, *local_spaces, out);
            } else {
                auto it = spaces.find(ownerrank);
                assert(it != spaces.end());
                std::shared_ptr<FESpacesAdapter> spaceptr = it->second;
                assert(std::distance(begin, end) > 0);
                write_element_selection(begin, end, *spaceptr, out);
            }

            CHECK_STREAM_WRITE_END("vol_proj", out);
        };

    };

    template class FESpaceSerializerDeserializer<2>;
    template class FESpaceSerializerDeserializer<3>;

    template<int Dimensions>
    class DefaultAlgorithm final : public TransferAssembler::Algorithm {
    public:
        using FunctionSpace = utopia::LibMeshFunctionSpace;
        using SparseMatrix  = utopia::USparseMatrix;
        using MeshBase      = libMesh::MeshBase;
        using DofMap        = libMesh::DofMap;
        using NTreeT 		= utopia::VTree<Dimensions>;
        using DataContainer = typename NTreeT::DataContainer;
        using Adapter       = typename NTreeT::DataType;

        DefaultAlgorithm(
            const std::shared_ptr<MeshBase> &from_mesh,
            const std::shared_ptr<DofMap>   &from_dofs,
            const std::shared_ptr<MeshBase> &to_mesh,
            const std::shared_ptr<DofMap>   &to_dofs,
            const TransferOptions &opts,
            const std::shared_ptr<LocalAssembler> &assembler,
            const std::shared_ptr<Local2Global>  &local2global)
        {
            init(from_mesh, from_dofs, to_mesh, to_dofs, opts, assembler, local2global);
        }

        void init(
            const std::shared_ptr<MeshBase> &from_mesh,
            const std::shared_ptr<DofMap>   &from_dofs,
            const std::shared_ptr<MeshBase> &to_mesh,
            const std::shared_ptr<DofMap>   &to_dofs,
            const TransferOptions &opts,
            const std::shared_ptr<LocalAssembler> &assembler,
            const std::shared_ptr<Local2Global>   &local2global)
        {
            this->from_mesh = from_mesh;
            this->from_dofs = from_dofs;
            this->to_mesh 	= to_mesh;
            this->to_dofs	= to_dofs;
            this->opts 		= opts;

            this->assembler = assembler;
            this->local2global = local2global;

            this->comm = moonolith::Communicator(from_mesh->comm().get());
            this->local_spaces = std::make_shared<FESpacesAdapter>(from_mesh, to_mesh, from_dofs, to_dofs, opts.from_var_num, opts.to_var_num);

            predicate = std::make_shared<moonolith::MasterAndSlave>();

            if(opts.tags.empty()){
                predicate->add(0, 1);
            } else {
                for(auto t : opts.tags) {
                    predicate->add(t.first, t.second);
                }
            }
        }

        inline bool assemble(Adapter &master,
                      Adapter &slave)
        {
            //FIXME assuming elements are all the same
             auto master_type = from_dofs->variable(opts.from_var_num).type();
             auto slave_type  = to_dofs->variable(opts.to_var_num).type();

            const auto &master_mesh = master.space();;
            const auto &slave_mesh  = slave.space();

            const int src_index  = master.element();
            const int dest_index = slave.element();

            auto &master_el = *master_mesh.elem(src_index);
            auto &slave_el  = *slave_mesh.elem(dest_index);

            return pg_assembler_.assemble(master_el, master_type, slave_el, slave_type,
                [&](std::vector<long> &master_dofs, std::vector<long> &slave_dofs) {
                    master_dofs = master.dof_map();
                    slave_dofs  = slave.dof_map();
                }
            );
        }

        bool assemble(std::vector<std::shared_ptr<SparseMatrix> > &mats)
        {
            if(assembler->n_forms() != mats.size()) {
                mats.resize(assembler->n_forms());
            }

            for(auto &mat_ptr : mats) {
                if(!mat_ptr) {
                    mat_ptr = std::make_shared<SparseMatrix>();
                }
            }

            mats_ = mats;
            return assemble_aux();
        }

        bool assemble_aux()
        {
            using namespace moonolith;

            init_tree();

            std::map<long, std::shared_ptr<FESpacesAdapter> > spaces;
            std::map<long, std::vector<std::shared_ptr<FESpacesAdapter> > > migrated_spaces;

            FESpaceSerializerDeserializer<Dimensions> serializer(
                from_mesh->comm(),
                opts,
                local_spaces);

            auto read = [&serializer] (
                const long ownerrank,
                const long senderrank,
                bool is_forwarding, DataContainer &data,
                InputStream &in
            ) {
                serializer.read(ownerrank, senderrank, is_forwarding, data, in);
            };

            auto write = [&serializer] (
                const long ownerrank, const long recvrank,
                const std::vector<long>::const_iterator &begin,
                const std::vector<long>::const_iterator &end,
                const DataContainer &data,
                OutputStream &out
            ) {
                serializer.write(ownerrank, recvrank, begin, end, data, out);
            };

            auto fun = [&](Adapter &master, Adapter &slave) -> bool {
                return this->assemble(master, slave);
            };

            pg_assembler_.initialize(
                comm,
                assembler,
                local2global,
                opts,
                from_dofs->n_dofs(),
                from_dofs->n_local_dofs(),
                to_dofs->n_dofs(),
                to_dofs->n_local_dofs()
            );

            moonolith::search_and_compute(comm, tree, predicate, read, write, fun, settings);

            pg_assembler_.finalize(mats_);
            pg_assembler_.print_stats();
            return true;
        }

        void init_tree()
        {
            using namespace moonolith;

            const auto n_elements_from = from_mesh->n_active_local_elem();
            const auto n_elements_to   = to_mesh->n_active_local_elem();
            const auto n_elements 	  = n_elements_from + n_elements_to;

            MOONOLITH_EVENT_BEGIN("create_adapters");

            tree = NTreeT::New(predicate, settings.max_elements, settings.max_depth);
            tree->reserve(n_elements);

            int offset = 0;
            if(opts.tags.empty()){
                int space_num = 0;

                for(auto s : local_spaces->spaces()) {
                    // const bool boundary_elements_only = opts.to_trace_space && space_num == 1;

                    if(s)
                    {
                        bool first = true;
                        libMesh::dof_id_type local_element_id = 0;
                        for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it, ++local_element_id) {

                            // if(boundary_elements_only) {
                            // 	if(!(*it)->is_boundary()) {
                            // 		continue;
                            // 	}
                            // }

                            auto elem = *it;
                            Adapter a(*s, elem->id(), offset+local_element_id,space_num);
                            assert(!local_spaces->dof_map(space_num)[local_element_id].empty());
                            a.set_dof_map(&local_spaces->dof_map(space_num)[local_element_id].global);
                            tree->insert(a);

                        }

                        offset += s->n_active_local_elem();
                    }

                    ++space_num;
                }

            } else {

                int space_num = 0;
                for(auto s : local_spaces->spaces()) {
                    if(s) {

                        bool first = true;

                        libMesh::dof_id_type local_element_id = 0;
                        for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it, ++local_element_id) {
                            auto elem=*it;
                            if (predicate->select(elem->subdomain_id())){
                                Adapter a(*s, elem->id(), offset+local_element_id,elem->subdomain_id());
                                assert(!local_spaces->dof_map(space_num)[local_element_id].empty());
                                a.set_dof_map(&local_spaces->dof_map(space_num)[local_element_id].global);
                                tree->insert(a);
                            }
                        }

                        offset += s->n_active_local_elem();
                    }

                    ++space_num;
                }
            }

            tree->root()->bound().static_bound().enlarge(1e-8);

            MOONOLITH_EVENT_END("create_adapters");
        }

    private:
        std::shared_ptr<MeshBase> from_mesh;
        std::shared_ptr<DofMap>   from_dofs;
        std::shared_ptr<MeshBase> to_mesh;
        std::shared_ptr<DofMap>   to_dofs;
        TransferOptions opts;

        std::shared_ptr<LocalAssembler> assembler;
        std::shared_ptr<Local2Global> local2global;

        moonolith::Communicator comm;
        moonolith::SearchSettings settings;
        std::shared_ptr<moonolith::MasterAndSlave> predicate;

        std::shared_ptr<NTreeT> tree;
        std::shared_ptr<FESpacesAdapter> local_spaces;

        private_::PetrovGalerkinAssembler pg_assembler_;
        std::vector<std::shared_ptr<SparseMatrix>> mats_;
    };

    template class DefaultAlgorithm<2>;
    template class DefaultAlgorithm<3>;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TransferAssembler::TransferAssembler(
        const std::shared_ptr<LocalAssembler> &assembler,
        const std::shared_ptr<Local2Global> &local2global)
    : assembler_(assembler), local2global_(local2global)
    {}

    TransferAssembler::~TransferAssembler() {}


    bool TransferAssembler::assemble(
        const std::shared_ptr<MeshBase> &from_mesh,
        const std::shared_ptr<DofMap>   &from_dofs,
        const std::shared_ptr<MeshBase> &to_mesh,
        const std::shared_ptr<DofMap>   &to_dofs,
        std::vector<std::shared_ptr<SparseMatrix> > &mats,
        const TransferOptions &opts)
    {
        assert(assembler_    && "assembler is required");
        assert(local2global_ && "local2global");

        ///////////////////////////

        moonolith::Communicator comm(from_mesh->comm().get());

        if(Utopia::instance().verbose()) {
            moonolith::root_describe("---------------------------------------\n"
                "begin: utopia::TransferAssembler::assemble",
                comm, std::cout);
        }

        Chrono c;
        c.start();

        ///////////////////////////
        // const int spatial_dimension = from_mesh->spatial_dimension();
        const int spatial_dimension = from_mesh->mesh_dimension();

        if(spatial_dimension == 1) {
            algorithm_ = std::make_shared<DefaultAlgorithm<1>>(from_mesh, from_dofs, to_mesh, to_dofs, opts, assembler_, local2global_);
        } else if(spatial_dimension == 2) {
            algorithm_ = std::make_shared<DefaultAlgorithm<2>>(from_mesh, from_dofs, to_mesh, to_dofs, opts, assembler_, local2global_);
        } else if(spatial_dimension == 3) {
            algorithm_ = std::make_shared<DefaultAlgorithm<3>>(from_mesh, from_dofs, to_mesh, to_dofs, opts, assembler_, local2global_);
        } else {
            assert(false && "dimension not supported");
            return false;
        }

        bool ok = algorithm_->assemble(mats);


        ///////////////////////////
        c.stop();

        if(Utopia::instance().verbose()) {
            std::stringstream ss;
            ss << "end: utopia::TransferAssembler::assemble\n";
            ss << c;
            ss << "---------------------------------------";
            ss << "\n";
            assembler_->print_stats(ss);
            moonolith::root_describe(ss.str(), comm, std::cout);
        }

        ///////////////////////////

        return ok;
    }

    bool TransferAssembler::assemble(
            const std::shared_ptr<MeshBase> &from_mesh,
            const std::shared_ptr<DofMap>   &from_dofs,
            const std::shared_ptr<MeshBase> &to_mesh,
            const std::shared_ptr<DofMap>   &to_dofs,
            SparseMatrix &B,
            const TransferOptions &opts
        )
    {
        std::vector<std::shared_ptr<SparseMatrix> > mats;
        mats.push_back(make_ref(B));
        return assemble(from_mesh, from_dofs, to_mesh, to_dofs, mats, opts);
    }
}

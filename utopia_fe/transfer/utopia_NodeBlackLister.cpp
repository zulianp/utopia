#include "utopia_NodeBlackLister.hpp"

namespace utopia {

    void NodeSetBlackLister::read(Input &in) 
    {
        in.get_all([this](Input &in) {
            int id = -1;

            in.get("id", id);

            if(id != -1) {
                node_sets_.insert(id);
            }
        });
    }

    bool NodeSetBlackLister::is_black_listed(const int node_set_id) const
    {
        return node_sets_.find(node_set_id) != node_sets_.end();
    }

    bool NodeSetBlackLister::is_black_listed(const libMesh::MeshBase &mesh, const libMesh::Elem &e) const
    {
        std::vector<libMesh::boundary_id_type> bids;
        const auto &info = mesh.get_boundary_info();

        std::size_t n_nodes = e.n_nodes();
        for(std::size_t i = 0; i < n_nodes; ++i) {
            const auto &n = e.node_ref(i);
            info.boundary_ids(&n, bids);

            for(auto id : bids) {
                if(is_black_listed(id)) {
                    return true;
                }
            }

        }

        return false;
    }

    bool NodeSetBlackLister::has_black_listed_nodes(const libMesh::MeshBase &mesh) const
    {
        const auto &bids = mesh.get_boundary_info().get_node_boundary_ids();

        for(auto i : bids) {
            if(is_black_listed(i)) {
                return true;
            }
        }

        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void ElementBlackList::read(Input &in) {
        black_lister_.read(in);
    }

    void ElementBlackList::init(const libMesh::MeshBase &mesh)
    {
        init(black_lister_, mesh);
    }

    bool ElementBlackList::is_black_listed(const SizeType element_id) {
        return black_listed_.find(element_id) != black_listed_.end();
    }

    bool ElementBlackList::is_black_listed(const libMesh::Elem &e) {
        if(is_boundary_mesh_) {
            auto parent = e.interior_parent();

            if(parent) {
                const bool found_black_listed = is_black_listed(parent->id());
                // if(found_black_listed) {
                //     std::cout << "FOUND BLACK LISTED" << std::endl;
                // }
                return found_black_listed;
            }
        } else {
            return is_black_listed(e.id());
        }

        return false;
    }

    ElementBlackList::ElementBlackList(const bool is_boundary_mesh)
    : is_boundary_mesh_(is_boundary_mesh)
    {}


    void ElementBlackList::init(const NodeSetBlackLister &black_lister, const libMesh::MeshBase &mesh)
    {
        if(black_lister.has_black_listed_nodes(mesh)) {
            for(auto eit = elements_begin(mesh); eit != elements_end(mesh); ++eit) {

                if(black_lister.is_black_listed(mesh, **eit)) {
                    black_listed_.insert((*eit)->id());
                }
            }
        }
    }

}

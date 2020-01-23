#ifndef UTOPIA_GRAPH_VIEW_HPP
#define UTOPIA_GRAPH_VIEW_HPP

#include "utopia_MeshView.hpp"
#include "utopia_StencilView.hpp"
#include "utopia_trilinos_Traits.hpp"

#include <memory>
#include <Tpetra_FECrsGraph.hpp>
#include <Tpetra_Map.hpp>

namespace utopia {

    template<class Mesh>
    class GraphView {};

    template<class Elem, class ExecutionSpace, typename...Args>
    class GraphView<TensorMeshView<Elem, ExecutionSpace, Args...>> {
    public:
        using NNZView             = Kokkos::View<std::size_t*, ExecutionSpace>;
        using OwnedPlusSharedView = Kokkos::View<SizeType*, ExecutionSpace>;


        NNZView nnz_;
    };


    template<class Mesh>
    class Graph {};


    template<class Elem, class Comm, class ExecutionSpace, typename...Args>
    class Graph<Mesh<Elem, Comm, ExecutionSpace, Uniform<Args...>>> {
    public:
        using Mesh        = utopia::Mesh<Elem, Comm, ExecutionSpace, Uniform<Args...>>;
        using MeshView    = typename Mesh::ViewDevice;
        using StencilView = utopia::StencilView<MeshView>;
        using Dev         = Traits<TpetraVector>::Device;
        using SizeType    = Traits<TpetraVector>::SizeType;
        using MapType     = TpetraVector::MapType;
        using GraphType   = Tpetra::FECrsGraph<SizeType, SizeType>;
        using View        = Kokkos::View<SizeType*, ExecutionSpace>;
        using DualNNZView = Kokkos::DualView<std::size_t*, ExecutionSpace>;

        void init(const Mesh &mesh)
        {
            auto mesh_view = mesh.view_device();
            StencilView stencil_view(mesh_view);

            std::cout <<  mesh.local_node_range() << std::endl;

            auto lnr = mesh.local_node_range();

            Dev::parallel_for(
                lnr,
                UTOPIA_LAMBDA(const SizeType &i) {

                std::cout << i << ") ";
                for(typename StencilView::SizeType k = 0; k < stencil_view.size(); ++k) {
                    auto idx = stencil_view.index(i, k);

                    if(idx != stencil_view.invalid()) {
                        std::cout << idx << " ";
                    }
                }

                std::cout << std::endl;

            });

            auto owned_row_map = Teuchos::rcp(
                new MapType(mesh.n_nodes(), lnr.extent(), 0, mesh.comm().get())
            );

            auto owned_plus_shared_row_map = owned_row_map;
            // auto owned_plus_shared_row_map = Teuchos::rcp(
            //     new MapType(
            //         Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
            //         index_list,
            //         0,
            //         mesh.comm().get()
            // ));

            auto graph =  Teuchos::rcp(new GraphType(
                owned_row_map,
                owned_plus_shared_row_map,
                stencil_view.size()
            ));

            // Range r(owned_row_map->getMinGlobalIndex(), owned_row_map->getMaxGlobalIndex() + 1);
            // SizeType cols[3];
            // for(SizeType i = r.begin(); i != r.end(); ++i) {
            //     cols[0] = i;

            //     if(i > 0) {
            //         cols[1] = i-1;
            //     } else {
            //         cols[1] = n_global-1;
            //     }

            //     if(i < n_global-1) {
            //         cols[2] = i + 1;
            //     } else {
            //         cols[2] = 0;
            //     }

            //     std::sort(std::begin(cols), std::end(cols));
            //     graph->insertGlobalIndices(i, 3, cols);
            // }

            // graph->fillComplete();
        }
    };

}

#endif //UTOPIA_GRAPH_VIEW_HPP

#ifndef UTOPIA_GRAPH_VIEW_HPP
#define UTOPIA_GRAPH_VIEW_HPP

#include "utopia_MeshView.hpp"
#include "utopia_StencilView.hpp"
#include "utopia_trilinos_Traits.hpp"

#include <memory>
#include <Tpetra_FECrsGraph.hpp>

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
        using MapType     = Tpetra::Map<SizeType, SizeType, ExecutionSpace>;
        using GraphType   = Tpetra::FECrsGraph<SizeType, SizeType>;
        using View        = Kokkos::View<SizeType*, ExecutionSpace>;
        using DualNNZView = Kokkos::DualView<std::size_t*, ExecutionSpace>;

        void init(const Mesh &mesh)
        {
            auto mesh_view = mesh.view_device();
            StencilView stencil_view(mesh_view);

            std::cout <<  mesh.local_node_range() << std::endl;

            Dev::parallel_for(
                mesh.local_node_range(),
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

            // Teuchos::RCP<const MapType> ownedRowMap, ownedPlusSharedRowMap;
            // // size_t maxNumEntriesPerRow = 3;

            // TpetraMatrix mat;

            // const SizeType rank = mat.comm().rank();
            // const SizeType size = mat.comm().size();

            // SizeType n_local = 10;
            // SizeType n_global = mat.comm().sum(n_local);

            // SizeType n_ghosts = (size>0) * 2;
            // View index_list("il", n_local + n_ghosts);
            // DualView nnz_z_row("nnz_z_row", n_local + n_ghosts);

            // auto nnz_z_row_dev = nnz_z_row.view_device();

            // Dev::parallel_for(n_local, UTOPIA_LAMBDA(const SizeType &i) {
            //     index_list(i) = i + n_local * rank;

            //     nnz_z_row_dev(i) = 3;

            //     if(i == 0 && size > 0) {
            //         index_list(n_local)     = rank == 0? n_global -1 : n_local * rank - 1;
            //         index_list(n_local + 1) = rank == size-1? 0 : n_local * (rank + 1);
            //         nnz_z_row_dev(n_local) = 3;
            //         nnz_z_row_dev(n_local + 1) = 3;
            //     }
            // });

            // ownedRowMap = Teuchos::rcp(new MapType(n_global, n_local, 0, mat.comm().get()));
            // ownedPlusSharedRowMap = Teuchos::rcp(
            //     new MapType(
            //         Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
            //         index_list,
            //         0,
            //         mat.comm().get()
            // ));

            // auto graph =  Teuchos::rcp(new GraphType(
            //     ownedRowMap,
            //     ownedPlusSharedRowMap,
            //     // maxNumEntriesPerRow
            //     nnz_z_row
            // ));

            // Range r(ownedRowMap->getMinGlobalIndex(), ownedRowMap->getMaxGlobalIndex() + 1);
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

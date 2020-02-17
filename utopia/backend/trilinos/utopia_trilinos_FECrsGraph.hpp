#ifndef UTOPIA_TRILINOS_FE_CRS_GRAPH_HPP
#define UTOPIA_TRILINOS_FE_CRS_GRAPH_HPP

#include "utopia_trilinos_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"
#include "utopia_trilinos_Traits.hpp"

#include <Tpetra_FECrsGraph.hpp>

namespace utopia {

    class TpetraFECrsGraph {
    public:
        //FIXME use correct execution space
        using Dev       = Traits<TpetraVector>::Device;
        using SizeType  = Traits<TpetraVector>::SizeType;
        using MapType   = Tpetra::Map<SizeType, SizeType>;
        using GraphType = Tpetra::FECrsGraph<SizeType, SizeType>;
        using View      = Kokkos::View<SizeType*>;
        using DualView  = Kokkos::DualView<std::size_t*>;

        TpetraFECrsGraph() {}

        inline const Teuchos::RCP<GraphType> &raw_type() const
        {
            return graph_;
        }

        inline Teuchos::RCP<GraphType> &raw_type()
        {
            return graph_;
        }

    private:
        Teuchos::RCP<GraphType> graph_;
    };

}

#endif //UTOPIA_TRILINOS_FE_CRS_GRAPH_HPP

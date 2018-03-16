#ifndef UTOPIA_TPETRAMATRIX_H
#define UTOPIA_TPETRAMATRIX_H

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"
#include "utopia_Size.hpp"

#include "utopia_Tpetra_Vector.hpp"

#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_Map_decl.hpp>

#include <iostream>

namespace utopia {
    
    class TpetraMatrix {
    public:
        
        typedef Tpetra::CrsMatrix<>                       crs_matrix_type;
        typedef Teuchos::RCP<crs_matrix_type>             rcp_crs_matrix_type;
        typedef Teuchos::RCP<const Teuchos::Comm<int> >   rcp_comm_type;
        typedef Tpetra::Map<>                             map_type;
        typedef Teuchos::RCP<map_type>                    rcp_map_type;
        typedef Tpetra::Vector<>::local_ordinal_type      local_ordinal_type;
        typedef Tpetra::Vector<>::global_ordinal_type     global_ordinal_type;
        typedef Tpetra::Vector<>::scalar_type             Scalar;
        
        TpetraMatrix() : owner_(true) {}
        
        // default serial constructor
        // TpetraMatrix (): owner_(true)
        // {
        //   auto map = Teuchos::rcp (new map_type ());
        //   mat_ = Teuchos::rcp (new crs_matrix_type (map, 0));
        // }
        
        //////////////////////////////////////////////////////////
        //constructor with only the communicator
        // TpetraMatrix(const rcp_comm_type & comm )
        // : _comm(comm), owner_(true)
        // {
        //   rcp_map_type _map = Teuchos::rcp (new map_type (_comm->getSize (), 0, _comm));
        //   rcp_crs_matrix_type mat_ (new crs_matrix_type (_map, 0));
        // }
        
        //////////////////////////////////////////////////////////
        // Map constructors based only on the given communicator//
        //////////////////////////////////////////////////////////
        
        //     template <class LocalOrdinal, class GlobalOrdinal, class Node>
        // TpetraMatrix (Tpetra::global_size_t numGlobalElem, GlobalOrdinal indexBase, const rcp_comm_type & comm, Tpetra::LocalGlobal lOrG, const Teuchos::RCP<Node> &node)
        // : _comm(comm), owner_(true)
        // {
        //   rcp_map_type _map (new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> (numGlobalElem, indexBase, &comm, lOrG, &node));
        //   rcp_crs_matrix_type mat_ (new crs_matrix_type (_map, 0));
        // }
        
        //     template <class LocalOrdinal, class GlobalOrdinal, class Node>
        // TpetraMatrix (Tpetra::global_size_t numGlobalElem, size_t numLocalElem, GlobalOrdinal indexBase, const rcp_comm_type &comm, const Teuchos::RCP<Node> &node)
        // : _comm(comm), owner_(true)
        // {
        //   rcp_map_type _map (new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> (numGlobalElem, numLocalElem, indexBase, &comm, &node));
        //   rcp_crs_matrix_type mat_ (new crs_matrix_type (_map, 0));
        // }
        
        //     template <class LocalOrdinal, class GlobalOrdinal, class Node>
        // TpetraMatrix (Tpetra::global_size_t numGlobalElem, const Teuchos::ArrayView<const GlobalOrdinal> &entryList, GlobalOrdinal indexBase, const rcp_comm_type &comm, const Teuchos::RCP<Node> &node)
        // : _comm(comm), owner_(true)
        // {
        //   rcp_map_type _map (new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> (numGlobalElem, &entryList, indexBase, &comm, &node));
        //   rcp_crs_matrix_type mat_ (new crs_matrix_type (_map, 0));
        // }
        
        //////////////////////////////////////////////////////////
        //constructor with a known map
        // TpetraMatrix(rcp_map_type  map)
        // {
        //   mat_.reset( new crs_matrix_type(map, 0));
        //   owner_ = true;
        // }
        
        /////////////////////////////////////////////////////////////
        // Matrix constructors based only on the given comminicator//
        /////////////////////////////////////////////////////////////
        
        //     template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
        //     class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
        //     class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
        //     class Node = Tpetra::Details::DefaultTypes::node_type,
        //     const bool classic = Node::classic>
        
        //         //RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>, DistObject<char, LocalOrdinal, GlobalOrdinal, Node, classic>
        //     TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
        //       size_t maxNumEntriesPerRow,
        //       Tpetra::ProfileType pftype = Tpetra::DynamicProfile,
        //       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
        //     : owner_(true)
        //     {
        //       mat_ = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (&rowMap, maxNumEntriesPerRow, pftype, &params));
        //     }
        // //
        //     template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
        //     class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
        //     class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
        //     class Node = Tpetra::Details::DefaultTypes::node_type,
        //     const bool classic = Node::classic>
        
        //     TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
        //       const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
        //       Tpetra::ProfileType pftype = Tpetra::DynamicProfile,
        //       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
        //     : owner_(true)
        //     {
        //       mat_ = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (&rowMap, &NumEntriesPerRowToAlloc, pftype, &params));
        //     }
        
        // //
        //     template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
        //     class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
        //     class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
        //     class Node = Tpetra::Details::DefaultTypes::node_type,
        //     const bool classic = Node::classic>
        
        //     TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
        //       const Teuchos::RCP<const map_type>& colMap,
        //       size_t maxNumEntriesPerRow,
        //       Tpetra::ProfileType pftype = Tpetra::DynamicProfile,
        //       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
        //     : owner_(true)
        //     {
        //       mat_ = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (&rowMap, &colMap, maxNumEntriesPerRow, pftype, &params));
        //     }
        
        
        // //
        //     template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
        //     class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
        //     class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
        //     class Node = Tpetra::Details::DefaultTypes::node_type,
        //     const bool classic = Node::classic>
        
        //     TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
        //       const Teuchos::RCP<const map_type>& colMap,
        //       const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
        //       Tpetra::ProfileType pftype = Tpetra::DynamicProfile,
        //       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
        //     : owner_(true)
        //     {
        //       mat_ = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
        //        (&rowMap, &colMap, &NumEntriesPerRowToAlloc, pftype, &params));
        
        //     }
        
        // //
        //     template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
        //     class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
        //     class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
        //     class Node = Tpetra::Details::DefaultTypes::node_type,
        //     const bool classic = Node::classic>
        
        //     TpetraMatrix (const Teuchos::RCP<const typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::crs_graph_type>& graph,
        //       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
        //     : owner_(true)
        //     {
        //       mat_ = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (&graph, &params));
        //     }
        
        // //
        //     template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
        //     class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
        //     class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
        //     class Node = Tpetra::Details::DefaultTypes::node_type,
        //     const bool classic = Node::classic>
        
        //     TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
        //       const Teuchos::RCP<const map_type>& colMap,
        //       const typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::row_map_type& rowPointers,
        //       const typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_graph_type::entries_type::non_const_type& columnIndices,
        //       const typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::values_type& values,
        //       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
        //     : owner_(true)
        //     {
        //       mat_ = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
        //        (&rowMap, &colMap, &rowPointers, &columnIndices, &values, &params));
        //     }
        
        // //
        //     template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
        //     class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
        //     class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
        //     class Node = Tpetra::Details::DefaultTypes::node_type,
        //     const bool classic = Node::classic>
        
        //     TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
        //       const Teuchos::RCP<const map_type>& colMap,
        //       const Teuchos::ArrayRCP<size_t>& rowPointers,
        //       const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
        //       const Teuchos::ArrayRCP<Scalar>& values,
        //       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
        //     : owner_(true)
        //     {
        //       mat_ = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
        //        (&rowMap, &colMap, &rowPointers, &columnIndices, &values, &params));
        //     }
        
        // //
        //     template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
        //     class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
        //     class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
        //     class Node = Tpetra::Details::DefaultTypes::node_type,
        //     const bool classic = Node::classic>
        
        //     TpetraMatrix(const Teuchos::RCP<const map_type>& rowMap,
        //      const Teuchos::RCP<const map_type>& colMap,
        //      const typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
        //      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
        //     : owner_(true)
        //     {
        //       mat_ = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
        //        (&rowMap, &colMap, &lclMatrix, &params));
        //     }
        
        /////////////////////////////////////////////////////////////
        ~TpetraMatrix()
        {}
        
        //deep copy
        //     template <class Node2>
        // rcp_crs_matrix_type  clone (
        //   const Teuchos::RCP<Node2> & node2,
        //   const Teuchos::RCP<Teuchos::ParameterList> & params = Teuchos::null)
        // {
        //   return mat_->clone(node2, params);
        // }
        
        // inline rcp_crs_matrix_type &implementation()
        // {
        //   return _mat;
        // }
        
        // inline const rcp_crs_matrix_type &implementation() const
        // {
        //   return _mat;
        // }
        
        // inline void set(global_ordinal_type gblRow,
        //  const Teuchos::ArrayView<global_ordinal_type> & cols,
        //  const Teuchos::ArrayView<scalar_type> & vals)
        // {
        //   mat_->insertGlobalValues(gblRow, cols, vals );
        // }
        
        inline void finalize()
        {
            mat_->fillComplete();
        }
        
        // inline global_ordinal_type global_index(local_ordinal_type row) const
        // {
        //     return mat_->getMap()->getGlobalElement(row);
        // }
        
        rcp_comm_type communicator() const
        {
            return mat_->getMap()->getComm();
        }
        
        void set_owner(const bool owner)
        {
            owner_ = owner;
        }
        
        
        //API functions
        
        void crs_init(const rcp_comm_type &comm,
                      std::size_t rows_local,
                      std::size_t cols_local,
                      Tpetra::global_size_t rows_global,
                      Tpetra::global_size_t cols_global,
                      std::size_t nnz_x_row);
        
        
        inline Range row_range() const
        {
            return  { mat_->getRowMap()->getMinGlobalIndex(), mat_->getRowMap()->getMaxGlobalIndex() + 1 };
        }
        
        inline Size size() const
        {
            return { mat_->getRowMap()->getGlobalNumElements(), mat_->getColMap()->getGlobalNumElements() };
        }
        
        inline Size local_size() const
        {
            return { mat_->getRowMap()->getNodeNumElements(), mat_->getColMap()->getNodeNumElements() };
        }
        
        inline void read_lock()
        {
            //TODO?
        }
        
        inline void read_unlock()
        {
            //TODO?
        }
        
        inline void write_lock()
        {
            //TODO?
        }
        
        inline void write_unlock()
        {
            this->finalize();
        }
        
        void describe(std::ostream &os = std::cout) const
        {
            auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(os));
            mat_->describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
        }
        
        inline void set(const global_ordinal_type &row, const global_ordinal_type &col, const Scalar &value)
        {
            mat_->insertGlobalValues(row, 1, &value, &col);
        }
        
        inline void add(const global_ordinal_type &row, const global_ordinal_type &col, const Scalar &value)
        {
            mat_->sumIntoGlobalValues(row, 1, &value, &col);
        }

        void mult(const TpetraVector &vec, TpetraVector &result) const;
        
    private:
        rcp_crs_matrix_type  mat_;
        bool                 owner_;
    }; //TpetraMatrix
}

#endif //UTOPIA_TPETRAMATRIX_H


#ifndef UTOPIA_TPETRAMATRIX_H
#define UTOPIA_TPETRAMATRIX_H

#include "utopia_trilinos.hpp"

#include "utopia_Base.hpp"

namespace utopia
{
typedef Tpetra::CrsMatrix<>                       crs_matrix_type;
typedef Teuchos::RCP<crs_matrix_type>             rcp_crs_matrix_type;
typedef Teuchos::RCP<const Teuchos::Comm<int> >   rcp_comm_type;
typedef Tpetra::Map<>                             map_type;
typedef Teuchos::RCP< const map_type >            rcp_map_type;
typedef Tpetra::Vector<>::local_ordinal_type      local_ordinal_type;
typedef Tpetra::Vector<>::global_ordinal_type     global_ordinal_type;
typedef Tpetra::Vector<>::scalar_type             scalar_type;

class TpetraMatrix
    {
    public:

        // default serial constructor
        TpetraMatrix (): _owner(true)
            {
            rcp_map_type _map = Teuchos::rcp (new map_type ());
            _comm = _map->getComm ();
            rcp_crs_matrix_type _mat (new crs_matrix_type (_map, 0));
            }

        //////////////////////////////////////////////////////////
        //constructor with only the communicator
        TpetraMatrix(const rcp_comm_type & comm )
            : _comm(comm), _owner(true)
            {
            rcp_map_type _map = Teuchos::rcp (new map_type (_comm->getSize (), 0, _comm));
            rcp_crs_matrix_type _mat (new crs_matrix_type (_map, 0));
            }

        //////////////////////////////////////////////////////////
        // Map constructors based only on the given communicator//
        //////////////////////////////////////////////////////////

        template <class LocalOrdinal, class GlobalOrdinal, class Node>
        TpetraMatrix (Tpetra::global_size_t numGlobalElem, GlobalOrdinal indexBase, const rcp_comm_type & comm, Tpetra::LocalGlobal lOrG, const Teuchos::RCP<Node> &node)
            : _comm(comm), _owner(true)
            {
            rcp_map_type _map (new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> (numGlobalElem, indexBase, &comm, lOrG, &node));
            rcp_crs_matrix_type _mat (new crs_matrix_type (_map, 0));
            }

        template <class LocalOrdinal, class GlobalOrdinal, class Node>
        TpetraMatrix (Tpetra::global_size_t numGlobalElements, size_t numLocalElem, GlobalOrdinal indexBase, const rcp_comm_type &comm, const Teuchos::RCP<Node> &node)
            : _comm(comm), _owner(true)
            {
            rcp_map_type _map (new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> (numGlobalElem, numLocalElements, indexBase, &comm, &node));
            rcp_crs_matrix_type _mat (new crs_matrix_type (_map, 0));
            }

        template <class LocalOrdinal, class GlobalOrdinal, class Node>
        TpetraMatrix (Tpetra::global_size_t numGlobalElem, const Teuchos::ArrayView<const GlobalOrdinal> &entryList, GlobalOrdinal indexBase, const rcp_comm_type &comm, const Teuchos::RCP<Node> &node)
            : _comm(comm), _owner(true)
            {
            rcp_map_type _map (new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> (numGlobalElem, &entryList, indexBase, &comm, &node));
            rcp_crs_matrix_type _mat (new crs_matrix_type (_map, 0));
            }

        //////////////////////////////////////////////////////////
        //constructor with a known map
        TpetraMatrix(rcp_map_type  map)
            {
            _comm = map->getComm ();
            _map = map;
            _mat.reset ( new crs_matrix_type (_map,0));
            _owner = true;
            }

        /////////////////////////////////////////////////////////////
        // Matrix constructors based only on the given comminicator//
        /////////////////////////////////////////////////////////////

        template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
                  class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
                  class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
                  class Node = Tpetra::Details::DefaultTypes::node_type,
                  const bool classic = Node::classic>

        //RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>, DistObject<char, LocalOrdinal, GlobalOrdinal, Node, classic>
        TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
                      size_t maxNumEntriesPerRow,
                      Tpetra::ProfileType pftype = Tpetra::DynamicProfile,
                      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
            : _owner(true)
            {
            rcp_crs_matrix_type _mat(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (&rowMap, maxNumEntriesPerRow, pftype, &params));
            _map = _mat->getMap();
            _comm = _map->getComm ();
            }
//
        template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
                  class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
                  class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
                  class Node = Tpetra::Details::DefaultTypes::node_type,
                  const bool classic = Node::classic>

        TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
                      const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
                      Tpetra::ProfileType pftype = Tpetra::DynamicProfile,
                      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
            : _owner(true)
            {
            rcp_crs_matrix_type _mat(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (&rowMap, &NumEntriesPerRowToAlloc, pftype, &params));
            _map = _mat->getMap();
            _comm = _map->getComm ();
            }

//
        template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
                  class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
                  class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
                  class Node = Tpetra::Details::DefaultTypes::node_type,
                  const bool classic = Node::classic>

        TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
                      const Teuchos::RCP<const map_type>& colMap,
                      size_t maxNumEntriesPerRow,
                      Tpetra::ProfileType pftype = Tpetra::DynamicProfile,
                      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
            : _owner(true)
            {
            rcp_crs_matrix_type _mat(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (&rowMap, &colMap, maxNumEntriesPerRow, pftype, &params));
            _map = _mat->getMap();
            _comm = _map->getComm ();
            }


//
        template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
                  class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
                  class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
                  class Node = Tpetra::Details::DefaultTypes::node_type,
                  const bool classic = Node::classic>

        TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
                      const Teuchos::RCP<const map_type>& colMap,
                      const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
                      Tpetra::ProfileType pftype = Tpetra::DynamicProfile,
                      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
            : _owner(true)
            {
            rcp_crs_matrix_type _mat(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
                                     (&rowMap, &colMap, &NumEntriesPerRowToAlloc, pftype, &params));
            _map = _mat->getMap();
            _comm = _map->getComm ();
            }

//
        template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
                  class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
                  class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
                  class Node = Tpetra::Details::DefaultTypes::node_type,
                  const bool classic = Node::classic>

        TpetraMatrix (const Teuchos::RCP<const typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::crs_graph_type>& graph,
                      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
            : _owner(true)
            {
            rcp_crs_matrix_type _mat(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (&graph, &params));
            _map = _mat->getMap();
            _comm = _map->getComm ();
            }

//
        template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
                  class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
                  class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
                  class Node = Tpetra::Details::DefaultTypes::node_type,
                  const bool classic = Node::classic>

        TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
                      const Teuchos::RCP<const map_type>& colMap,
                      const typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::row_map_type& rowPointers,
                      const typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_graph_type::entries_type::non_const_type& columnIndices,
                      const typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::values_type& values,
                      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
            : _owner(true)
            {
            rcp_crs_matrix_type _mat(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
                                     (&rowMap, &colMap, &rowPointers, &columnIndices, &values, &params));
            _map = _mat->getMap();
            _comm = _map->getComm ();
            }

//
        template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
                  class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
                  class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
                  class Node = Tpetra::Details::DefaultTypes::node_type,
                  const bool classic = Node::classic>

        TpetraMatrix (const Teuchos::RCP<const map_type>& rowMap,
                      const Teuchos::RCP<const map_type>& colMap,
                      const Teuchos::ArrayRCP<size_t>& rowPointers,
                      const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
                      const Teuchos::ArrayRCP<Scalar>& values,
                      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
            : _owner(true)
            {
            rcp_crs_matrix_type _mat(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
                                     (&rowMap, &colMap, &rowPointers, &columnIndices, &values, &params));
            _map = _mat->getMap();
            _comm = _map->getComm ();
            }

//
        template <class Scalar = Tpetra::Details::DefaultTypes::scalar_type,
                  class LocalOrdinal = Tpetra::Details::DefaultTypes::local_ordinal_type,
                  class GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
                  class Node = Tpetra::Details::DefaultTypes::node_type,
                  const bool classic = Node::classic>

        TpetraMatrix(const Teuchos::RCP<const map_type>& rowMap,
                     const Teuchos::RCP<const map_type>& colMap,
                     const typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
                     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
            : _owner(true)
            {
            rcp_crs_matrix_type _mat(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
                                     (&rowMap, &colMap, &lclMatrix, &params));
            _map = _mat->getMap();
            _comm = _map->getComm ();
            }

        /////////////////////////////////////////////////////////////
        ~TpetraMatrix()
            {
            if(_owner)
                {
                _mat->~CrsMatrix();
                }
            }

        //deep copy
        template <class Node2>
        rcp_crs_matrix_type  clone (
            const Teuchos::RCP<Node2> & node2,
            const Teuchos::RCP<Teuchos::ParameterList> & params = Teuchos::null)
            {
            return _mat->clone(node2, params);
            }

        inline rcp_crs_matrix_type &implementation()
            {
            return _mat;
            }

        inline const rcp_crs_matrix_type &implementation() const
            {
            return _mat;
            }

        inline bool insertGlobalValues(global_ordinal_type gblRow,
                                       const Teuchos::ArrayView<global_ordinal_type> & cols,
                                       const Teuchos::ArrayView<scalar_type> & vals)
            {
            _mat->insertGlobalValues(gblRow, cols, vals );
            return true;
            }

        inline bool fillComplete()
            {
            _mat->fillComplete();
            return true;
            }

        inline size_t getNodeNumElements()
            {
            return _map->getNodeNumElements();
            }

        inline bool getGlobalElement(local_ordinal_type row)
            {
            _map->getGlobalElement(row);
            return true;
            }

        rcp_comm_type &communicator()
            {
            return _comm;
            }

        void set_owner(const bool owner)
            {
            _owner = owner;
            }

    protected:
        rcp_crs_matrix_type &operator= (const rcp_crs_matrix_type &matrix) { }

    private:
        rcp_comm_type        _comm;
        rcp_map_type         _map;
        rcp_crs_matrix_type  _mat;
        bool                 _owner;
    }; //TpetraMatrix
}

#endif //UTOPIA_TPETRAMATRIX_H

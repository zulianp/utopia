#ifndef UTOPIA_STK_FORWARD_DECLARATIONS_HPP
#define UTOPIA_STK_FORWARD_DECLARATIONS_HPP

namespace stk {

    namespace mesh {

        class BulkData;
        class MetaData;
        class Selector;
        class Bucket;

    }  // namespace mesh

    namespace io {
        class StkMeshIoBroker;
    }
}  // namespace stk

namespace utopia {

    namespace stk {

        class Mesh;
        class FunctionSpace;
        class DofMap;
        class MeshIO;
        class SpaceIO;
        // class StkIntrepid2Assembler;
        class Transport;
        class Mass;

    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_FORWARD_DECLARATIONS_HPP
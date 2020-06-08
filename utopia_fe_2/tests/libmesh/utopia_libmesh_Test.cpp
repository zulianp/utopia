#include "utopia.hpp"
// #include "utopia_ParallelTestRunner.hpp"
#include "utopia_Testing.hpp"

#include "utopia_libmesh.hpp"

namespace utopia {

    // FIXME
    template <class Mesh>
    class MeshTest final /*: public UnitTest<PetscCommunicator> */ {
    public:
        using Comm = typename Traits<Mesh>::Communicator;
        void run() /*override*/ {
            Comm comm;

            InputParameters params;

            params.set("format", "libmesh");
            params.set("mesh_type", "square");
            params.set("elem_type", "TRI3");

            Mesh mesh(comm);
            mesh.read(params);

            mesh.describe(std::cout);
        }
    };

    template <class Mesh>
    class FunctionSpaceTest final /*: public UnitTest<PetscCommunicator> */ {
    public:
        using Comm = typename Traits<Mesh>::Communicator;
        void run() /*override*/ {
            Comm comm;

            InputParameters params;

            params.set("format", "libmesh");
            params.set("mesh_type", "cube");
            params.set("elem_type", "TET4");

            FunctionSpace<Mesh> space(comm);
            space.read(params);

            space.describe(std::cout);
        }
    };

    static void libmesh_specific() {
        // FIXME
        MeshTest<LMMesh>().run();
        FunctionSpaceTest<LMMesh>().run();
    }

    UTOPIA_REGISTER_TEST_FUNCTION(libmesh_specific);
}  // namespace utopia

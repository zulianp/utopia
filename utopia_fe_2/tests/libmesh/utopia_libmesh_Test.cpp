#include "utopia.hpp"
// #include "utopia_ParallelTestRunner.hpp"
#include "utopia_Testing.hpp"

#include "utopia_libmesh.hpp"

namespace utopia {

    // FIXME
    template <class Mesh>
    class MeshTest final /*: public UnitTest<PetscCommunicator> */ {
    public:
        void run() /*override*/ { disp("HI: MeshTest"); }
    };

    static void libmesh_specific() {
        // FIXME
        MeshTest<LMMesh>().run();
        // run_parallel_test<MeshTest<LMMesh>>();
    }

    UTOPIA_REGISTER_TEST_FUNCTION(libmesh_specific);
}  // namespace utopia

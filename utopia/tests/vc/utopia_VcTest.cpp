#include "utopia_Testing.hpp"
#include "utopia_Vc.hpp"

namespace utopia {

    class VcTest {
    public:
        void run() {}
    };

    static void Vc_ops() { VcTest().run(); }

    UTOPIA_REGISTER_TEST_FUNCTION(Vc_ops);
}  // namespace utopia
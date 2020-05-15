#ifndef UTOPIA_INTREPID2_TEST_HPP
#define UTOPIA_INTREPID2_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace libMesh {
    class LibMeshInit;
}

namespace utopia {
    //  void run_intrepid2_test(libMesh::LibMeshInit &init);

    class Intrepid2Test final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "intrepid2"; }
    };

}  // namespace utopia

#endif  // UTOPIA_INTREPID2_TEST_HPP

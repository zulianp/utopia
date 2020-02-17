#ifndef UTOPIA_INTREPID2_TEST_HPP
#define UTOPIA_INTREPID2_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace libMesh {
    class LibMeshInit;
}


namespace utopia {
  //  void run_intrepid2_test(libMesh::LibMeshInit &init);

    class Intrepid2Test final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "intrepid2";
        }
    };

}



#endif //UTOPIA_INTREPID2_TEST_HPP

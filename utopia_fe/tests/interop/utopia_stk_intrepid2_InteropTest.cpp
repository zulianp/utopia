#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"

#include "utopia_FEInteroperability.hpp"

#include "utopia_intrepid2_LaplaceOperator.hpp"

using namespace utopia;

void poisson_problem() {
    using FunctionSpace_t = utopia::stk::FunctionSpace;
    using FE_t = utopia::intrepid2::FE<double>;

    auto params = param_list(param("path", "../data/knf/pump/membrane.e"));

    FunctionSpace_t space;
    space.read(params);

    // Chrono c;
    // c.start();

    auto fe_ptr = std::make_shared<FE_t>();
    create_fe(space, *fe_ptr, 0);

    LaplaceOperator<double> lapl{1.0};
    intrepid2::Assemble<LaplaceOperator<double>> assembler(lapl, fe_ptr);
    assembler.init();

    // c.stop();
    // std::cout << c << std::endl;

    // std::ofstream os("prova.txt");
    // assembler.describe(os);
    // os.close();
}

void interop_stk_intrepid2() { UTOPIA_RUN_TEST(poisson_problem); }

UTOPIA_REGISTER_TEST_FUNCTION(interop_stk_intrepid2);

#endif  // UTOPIA_WITH_INTREPID2
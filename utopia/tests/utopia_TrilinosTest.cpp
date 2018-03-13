#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS

#include "utopia_TrilinosTest.hpp"
#include "utopia.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_solvers.hpp"

namespace utopia {   

    void trilinos_build_test()
    {   
        auto n = 10;
        TVectord v = local_values(n, 1.);

        //FIXME replace this with an actual test
        disp(v);
    }

    void trilinos_axpy_test()
    {
       auto n = 10;
       TVectord y = local_values(n, 1.);
       TVectord x = local_values(n, 5.);
       auto alpha = 0.1; 
       y += alpha * x;

       //FIXME replace this with an actual test
       disp(y);
    }

    void run_trilinos_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("TrilinosTest");
        UTOPIA_RUN_TEST(trilinos_build_test);
        UTOPIA_RUN_TEST(trilinos_axpy_test);
        UTOPIA_UNIT_TEST_END("TrilinosTest");
    }
}

#else //WITH_TRILINOS

namespace utopia
{
    void run_trilinos_test() {}
}

#endif //WITH_TRILINOS

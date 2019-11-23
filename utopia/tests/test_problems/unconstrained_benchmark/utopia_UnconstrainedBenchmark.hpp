#ifndef UTOPIA_UNCONSTRAINED_TEST_SET_INCLUDES
#define UTOPIA_UNCONSTRAINED_TEST_SET_INCLUDES

    #include "utopia_TestFunctions.hpp"

    // Ordering of function follows paper:
    // More, Garbow, Hillstrom - Testing unconstrained optimization software - Unconstrained minimization section
    #include "utopia_01Rosenbrock.hpp"
    #include "utopia_03Powell.hpp"
    #include "utopia_04Brown.hpp"
    #include "utopia_05Beale.hpp"
    #include "utopia_07Helical.hpp"
    #include "utopia_09Gaussian.hpp"
    #include "utopia_11Gulf.hpp"
    #include "utopia_12Box.hpp"
    #include "utopia_14Woods.hpp"
    #include "utopia_16BrownDennis.hpp"
    #include "utopia_18Biggs.hpp"
    #include "utopia_20Watson.hpp"
    #include "utopia_21ExtendedRosenbrock.hpp"
    #include "utopia_22ExtendedPowellSingular.hpp"
    #include "utopia_23Penalty1.hpp"
    #include "utopia_24Penalty2.hpp"
    #include "utopia_25VariablyDim.hpp"
    #include "utopia_26Trigonometric.hpp"
    #include "utopia_35Chebyquad.hpp"
    
    // Some other test functions
    // TODO::make them follow interface 
    #include "utopia_QPTestFunction2D.hpp"
    #include "utopia_RastriginTestFunction.hpp"
    #include "utopia_TestFunctionsND.hpp"

#endif //UTOPIA_UNCONSTRAINED_TEST_SET_INCLUDES

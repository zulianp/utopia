#ifndef UTOPIA_DESCRIBE_FUNCTION_SPACE_APP
#define UTOPIA_DESCRIBE_FUNCTION_SPACE_APP

#include "utopia_fe_Core.hpp"
#include "utopia_fe_base.hpp"

#include "utopia_Input.hpp"

namespace utopia {

    template <class FunctionSpace>
    class DescribeFunctionSpaceApp {
    public:
        static void run(Input &in) {
            FunctionSpace space;
            space.read(in);

            if (space.empty()) {
                space.comm().root_print("No input space and mesh provided, using unit cube!\n");

                space.mesh().unit_cube(3, 3, 3);
                space.initialize();
            }

            space.describe(utopia::out().stream());
        }
    };

}  // namespace utopia

#endif  // UTOPIA_DESCRIBE_FUNCTION_SPACE_APP


#include "utopia_Testing.hpp"

#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"
#include "utopia_ui.hpp"

#include "utopia_RunParallelTest.hpp"

#include "utopia_Traits.hpp"

#include "utopia_Testing.hpp"
#include "utopia_UnitTest.hpp"

namespace utopia {

    template <class FunctionSpace>
    class FunctionSpaceTest final : public UnitTest<typename Traits<FunctionSpace>::Communicator> {
    public:
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using SizeType = typename Traits<FunctionSpace>::SizeType;
        using Vector = typename Traits<FunctionSpace>::Vector;

        void declare_field() {
            auto params = param_list(param("mesh", param_list(param("type", "cube"))));
            FunctionSpace space;
            space.read(params);
            space.template declare_new_nodal_field<Scalar>("velocity", 3);

            Vector v, v2;
            space.create_vector(v);
            space.create_vector(v2);

            v.set(1.0);
            space.global_vector_to_nodal_field(v);
            space.nodal_field_to_global_vector(v2);
        }

        void run() override { UTOPIA_RUN_TEST(declare_field); }
    };

}  // namespace utopia

using namespace utopia;

void ustk_fs() { utopia::run_parallel_test<FunctionSpaceTest<utopia::stk::FunctionSpace>>(); }

UTOPIA_REGISTER_TEST_FUNCTION(ustk_fs);

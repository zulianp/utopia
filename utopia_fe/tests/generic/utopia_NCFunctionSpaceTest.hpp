#ifndef UTOPIA_NC_FUNCTION_SPACE_TEST_HPP
#define UTOPIA_NC_FUNCTION_SPACE_TEST_HPP

#include "utopia_Testing.hpp"
#include "utopia_Traits.hpp"
#include "utopia_UnitTest.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_NCFunctionSpace.hpp"

namespace utopia {

    template <class FunctionSpace>
    class NCFunctionSpaceTest final : public UnitTest<typename Traits<FunctionSpace>::Communicator> {
    public:
        using SizeType = typename Traits<FunctionSpace>::SizeType;
        using Comm = typename Traits<FunctionSpace>::Communicator;
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Mesh = typename Traits<FunctionSpace>::Mesh;
        using NCFunctionSpace = utopia::NCFunctionSpace<FunctionSpace>;

#ifdef UTOPIA_ENABLE_MOONOLITH
        void mortar() {
            // FIXME
            if(Comm::get_default().size() > 1) return;

            NCFunctionSpace ncspace(this->comm());
            ncspace.import("../data/testing/decomposition.yaml");

            Vector v, vi;
            ncspace.create_vector(v);
            v.set(1.);

            ncspace.project(v, vi);

            ncspace.write("interp.e", vi);
        }
#endif

        void run() override {
#ifdef UTOPIA_ENABLE_MOONOLITH
            UTOPIA_RUN_TEST(mortar);
#endif  // UTOPIA_ENABLE_MOONOLITH
        }
    };

}  // namespace utopia

#endif  // UTOPIA_NC_FUNCTION_SPACE_TEST_HPP
#ifndef UTOPIA_SCALAR_PRODUCT_TEST_HPP
#define UTOPIA_SCALAR_PRODUCT_TEST_HPP

#include "utopia_Testing.hpp"
#include "utopia_UnitTest.hpp"

#include "utopia_fe_base.hpp"
#include "utopia_ui.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_intrepid2_L2Norm.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ScalarProductTest {
    public:
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using SideSet_t = typename Traits<FunctionSpace>::SideSet;
        using Field_t = utopia::Field<FunctionSpace>;

        InputParameters cube_space_param(const int n_var) const {
            const int nx = 2;
            const int ny = 2;
            const int nz = 2;

            return param_list(
                param("n_var", n_var),
                param("mesh", param_list(param("type", "cube"), param("nx", nx), param("ny", ny), param("nz", nz))));
        }

        // void uniform_function_on_cube() {
        //     FunctionSpace space;

        //     int n_var = 3;
        //     auto params = cube_space_param(n_var);
        //     space.read(params);

        //     Field<FunctionSpace> field;
        //     space.create_field(field);

        //     field.data().set(1.0);

        //     std::vector<Scalar_t> norms;
        //     intrepid2::l2_norm(field, norms);

        //     UTOPIA_TEST_EQ(n_var, norms.size());

        //     for (int c = 0; c < n_var; ++c) {
        //         utopia_test_assert(approxeq(1.0, norms[c], 1e-8));
        //     }
        // }

        void run() {
            // UTOPIA_RUN_TEST(uniform_function_on_cube);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_SCALAR_PRODUCT_TEST_HPP

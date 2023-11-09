
#include "utopia_Testing.hpp"

#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"
#include "utopia_stk_SpaceIO.hpp"
#include "utopia_ui.hpp"

#include "utopia_RunParallelTest.hpp"

#include "utopia_Traits.hpp"

#include "utopia_Testing.hpp"
#include "utopia_UnitTest.hpp"

namespace utopia {

    template <class FunctionSpace>
    class RestartTest final : public UnitTest<typename Traits<FunctionSpace>::Communicator> {
    public:
        using Mesh = typename Traits<FunctionSpace>::Mesh;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using SizeType = typename Traits<FunctionSpace>::SizeType;
        using Vector = typename Traits<FunctionSpace>::Vector;

        void test_write_and_read_back() {
            auto mesh = std::make_shared<Mesh>();
            mesh->unit_cube(6, 6, 6);

            FunctionSpace space(mesh);
            space.initialize();

            Field<FunctionSpace> x("x"), y("y"), z("z"), p("p");
            space.create_nodal_vector_field(3, x);
            space.create_nodal_vector_field(3, y);
            space.create_nodal_vector_field(3, z);
            space.create_nodal_vector_field(1, p);

            double t = 1;

            {  // Write
                IO<FunctionSpace> output(space);
                output.set_output_path("./restart_test.e");

                // To be done only once!
                output.register_output_field(x);
                output.register_output_field(y);
                output.register_output_field(z);
                output.register_output_field(p);

                x.data().set(0.0);
                y.data().set(1.0);
                z.data().set(2.0);
                p.data().set(3.0);

                // To be done every timestep
                output.update_output_field(x);
                output.update_output_field(y);
                output.update_output_field(z);
                output.update_output_field(p);
                utopia_test_assert(output.write(0, t));
            }

            {  // Read
                IO<FunctionSpace> input(space);
                input.set_input_path("./restart_test.e");
                input.import_all_field_data(true);
                input.open_input();
                input.load_time_step(t);

                Field<FunctionSpace> x_disk("x", make_ref(space)), y_disk("y", make_ref(space)),
                    z_disk("z", make_ref(space)), p_disk("p", make_ref(space));

                utopia_test_assert(input.read_nodal(x_disk));
                utopia_test_assert(input.read_nodal(y_disk));
                utopia_test_assert(input.read_nodal(z_disk));
                utopia_test_assert(input.read_nodal(p_disk));

                Scalar diff_x = norm2(x_disk.data() - x.data());
                Scalar diff_y = norm2(y_disk.data() - y.data());
                Scalar diff_z = norm2(z_disk.data() - z.data());
                Scalar diff_p = norm2(p_disk.data() - p.data());

                // test that fields are the same!
                utopia_test_assert(diff_x == 0);
                utopia_test_assert(diff_y == 0);
                utopia_test_assert(diff_z == 0);
                utopia_test_assert(diff_p == 0);
            }
        }

        void run() override { UTOPIA_RUN_TEST(test_write_and_read_back); }
    };

}  // namespace utopia

using namespace utopia;

void stk_restart() { utopia::run_parallel_test<RestartTest<utopia::stk::FunctionSpace>>(); }

UTOPIA_REGISTER_TEST_FUNCTION(stk_restart);

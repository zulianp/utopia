#include "utopia.hpp"
#include "utopia_IsSubTree.hpp"
#include "utopia_Testing.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#ifdef UTOPIA_ENABLE_TRILINOS
#include "utopia_trilinos.hpp"
#endif

namespace utopia {

    template <class Matrix, class Vector>
    class DeviceTest {
    public:
        using Traits = utopia::Traits<Vector>;
        using SizeType = typename Traits::SizeType;
        using Scalar = typename Traits::Scalar;
        using Comm = typename Traits::Communicator;

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << Traits::backend_info().get_name() << std::endl;
            }
        }

        static void test_range_device() {
            Comm world;
            // host context
            const SizeType n_local = 2;
            const SizeType n = n_local * world.size();

            Vector v(layout(world, n_local, n), 1.);

            {
                // device context
                auto r = range_device(v);
                auto v_device = view_device(v);

                // using kokkos parallel for (r contains the traits of the backend
                // and it should be used to provide specific range policies associated to utopia types
                parallel_for(
                    r, UTOPIA_LAMBDA(const SizeType &i) { v_device.set(i, i); });

                Scalar reduce_v = 0.0;
                parallel_reduce(
                    r, UTOPIA_LAMBDA(const SizeType &i) { return v_device.get(i); }, reduce_v);
            }

            // host context
            SizeType sum_v = sum(v);
            utopia_test_assert(sum_v == ((n - 1) * (n)) / 2);
        }

        static void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(test_range_device);
        }
    };

    static void device_test() {
#ifdef UTOPIA_ENABLE_BLAS
        DeviceTest<BlasMatrixd, BlasVectord>().run();
#endif  // UTOPIA_ENABLE_BLAS

#ifdef UTOPIA_ENABLE_PETSC
        DeviceTest<utopia::PetscMatrix, utopia::PetscVector>::run();
#endif  // UTOPIA_ENABLE_PETSC

#ifdef UTOPIA_ENABLE_TRILINOS
        DeviceTest<utopia::TpetraMatrix, utopia::TpetraVector>::run();
#endif  // UTOPIA_ENABLE_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(device_test);

}  // namespace utopia
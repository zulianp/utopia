#include "utopia.hpp"
#include "utopia_IsSubTree.hpp"
#include "utopia_Testing.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#ifdef WITH_TRILINOS
#include "utopia_trilinos.hpp"
#endif

namespace utopia {

    template<class Matrix, class Vector>
    class DeviceTest {
    public:
        using Traits   = utopia::Traits<Vector>;
        using SizeType = typename Traits::SizeType;
        using Scalar   = typename Traits::Scalar;

        static void print_backend_info()
        {
            if(Utopia::instance().verbose() && mpi_world_rank() == 0) {
                std::cout << "\nBackend: " << Traits::backend_info().get_name() << std::endl;
            }
        }

        static void test_range_device()
        {
            //host context
            SizeType n = 10;
            Vector v = values(n, 1.);

            {
                //device context
                auto r        = v.range_device();
                auto v_device = view_device(v);

                //using kokkos parallel for (r contains the traits of the backend
                //and it should be used to provide specific range policies associated to utopia types
                parallel_for(r, UTOPIA_LAMBDA(const SizeType &i) {
                    v_device.set(i, i);
                });

                Scalar reduce_v = 0.0;
                parallel_reduce(r, UTOPIA_LAMBDA(const SizeType &i) {
                    return v_device.get(i);
                }, reduce_v);


                //test local reduce
                reduce_v = v.comm().sum(reduce_v);
                utopia_test_assert(static_cast<SizeType>(reduce_v) == (n - 1) * (n/2));
            }

            //host context
            SizeType sum_v = sum(v);
            utopia_test_assert(sum_v == (n - 1) * (n/2));
        }

        static void run()
        {
            print_backend_info();
            UTOPIA_RUN_TEST(test_range_device);
        }
    };

    static void device_test()
    {
// #ifdef WITH_BLAS
//         DeviceTest<BlasMatrixd, BlasVectord>().run();
// #endif //WITH_BLAS

#ifdef WITH_PETSC
        DeviceTest<PetscMatrix, PetscVector>().run();
#endif //WITH_PETSC

#ifdef WITH_TRILINOS
        DeviceTest<TpetraMatrix, TpetraVector>().run();
#endif //WITH_TRILINOS

    }

    UTOPIA_REGISTER_TEST_FUNCTION(device_test);

}
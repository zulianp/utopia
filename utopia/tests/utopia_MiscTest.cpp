#include "utopia.hpp"
#include "utopia_MiscTest.hpp"
#include "utopia_FactoryMethod.hpp"

namespace utopia {

    #ifdef WITH_CUDA
    void cuda_hello_world()
    {
        //For the moment no hybrid
        if(utopia::mpi_world_size() != 1) return;

        const int n = 100;
        CUMatrixd m = identity(n, n);
        CUVectord v = values(n, 1.0/3.0);
        const double r = dot(v, v);
        CUVectord mv = m * v;

        const double mvr =  dot(mv, mv);

        // disp(r);
        // disp(mvr);
    }
    #endif


#ifdef WITH_UTOPIA_OPENCL
    void test_opencl_code()
    {
        const int n = 30;

        CLMatrixd m  = identity(n, n);
        CLVectord v1 = values(n, 0.1),
                  v2 = values(n, 0.2),
                  v3 = values(n, 0.3),
                  v4 = values(n, 0.4);

        CLVectord res = abs(m * pow2(v1) + sqrt(v2) - v3);

        // std::cout << "-------------------\n";
        // disp(v1);
        // std::cout << "-------------------\n";
        // disp(res);
        // std::cout << "-------------------\n";

        CLMatrixd mat_res = transpose(abs(0.1 * (m * m) - m) * (m)); //inner transpose does not work yet.

        // std::cout << "-------------------\n";
        // // disp(mat_res);
        // std::cout << "-------------------\n";
    }

#endif //WITH_UTOPIA_OPENCL

#ifdef WITH_LAPACK
#ifdef WITH_BLAS
    void test_lapack_eigen_solver()
    {
        using namespace utopia;
        using SizeType = Traits<BlasMatrixd>::SizeType;
        
        const int n = 9;
        BlasMatrixd A = zeros(n, n);
        std::vector<SizeType> is{0, 1, 2, 3, 4, 5, 6, 7, 8};

        {
            Write<BlasMatrixd> w_A(A);
            A.set_matrix(is, is, {
                +6.944444e+05, 0, 0, 0, 0, 0, 0, 0, 0,
                0, +1.388889e+06, -5.555556e+04, 0, -3.194444e+05, -3.472222e+05, 0, 0, 0,
                0, -5.555556e+04, +1.388890e+06, 0, -3.472222e+05, -3.194448e+05, -5.555555e-02, 0, -3.472222e-01,
                0, 0, 0, +6.944451e+05, 0, 0, 0, 0, 0,
                0, -3.194444e+05, -3.472222e+05, 0, +6.944444e+05, -2.777778e+04, 0, 0, 0,
                0, -3.472222e+05, -3.194448e+05, 0, -2.777778e+04, +6.944451e+05, -3.472222e-01, 0, -2.777778e-02,
                0, 0, -5.555555e-02, 0, 0, -3.472222e-01, +1.388889e+00, 0, -3.194444e-01,
                0, 0, 0, 0, 0, 0, 0, +6.944444e-01, 0,
                0, 0, -3.472222e-01, 0, 0, -2.777778e-02, -3.194444e-01, 0, +6.944444e-01
            });
        }

        BlasMatrixd B = zeros(n, n);
        {
            Write<BlasMatrixd> w_B(B);
            B.set_matrix(is, is, {
                +6.944444e+05, 0, 0, 0, 0, 0, 0, 0, 0,
                0, +2.387253e+06, 0, 0, 0, 0, 0, 0, 0,
                0, 0, +2.387802e+06, 0, 0, 0, 0, 0, 0,
                0, 0, 0, +6.944451e+05, 0, 0, 0, 0, 0,
                0, 0, 0, 0, +1.193627e+06, 0, 0, 0, 0,
                0, 0, 0, 0, 0, +1.193901e+06, 0, 0, 0,
                0, 0, 0, 0, 0, 0, +1.841198e+00, 0, 0,
                0, 0, 0, 0, 0, 0, 0, +6.944444e-01, 0,
                0, 0, 0, 0, 0, 0, 0, 0, +9.205991e-01
            });

         }

        BlasMatrixd V;

        {
            Write<BlasMatrixd> w(A);
            A.set(0, 1, 2);
            A.set(1, 0, 2);
        }

        BlasVectord e;

        bool ok;
        ok = spd_geig_small(A, B, 0.3, e, V); utopia_test_assert(ok);
        ok = spd_geig(A, B, e, V);            utopia_test_assert(ok);
        ok = spd_eig(A, e, V);                utopia_test_assert(ok);
    }

    void factory_method()
    {
        FactoryMethod<BlasMatrixd> factory;
        auto ptr = factory.make();

        auto fact_1 = make_factory<BlasMatrixd, BlasMatrixd>();
        auto ptr_2  = fact_1->make();
    }

#endif //WITH_BLAS
#endif //WITH_LAPACK


    void runMiscTest() {
        // std::cout << "Begin: MiscTest" << std::endl;
        UTOPIA_UNIT_TEST_BEGIN("MiscTest");

#ifdef WITH_CUDA
           UTOPIA_RUN_TEST(cuda_hello_world);
#endif

#ifdef WITH_UTOPIA_OPENCL
            UTOPIA_RUN_TEST(test_opencl_code);
#endif

#ifdef WITH_LAPACK
#ifdef WITH_BLAS
            UTOPIA_RUN_TEST(test_lapack_eigen_solver);
#endif
#endif
        UTOPIA_UNIT_TEST_END("MiscTest");

    }
}

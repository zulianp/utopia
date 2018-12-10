#include "utopia.hpp"
#include "utopia_AlgebraTest.hpp"
#include "utopia_IsSubTree.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class AlgebraTest {
    private:
        typedef typename utopia::Traits<Vector>::Scalar Scalar;

        void norm_test()
        {
            Vector v = zeros(2);
            {
                auto r = range(v);
                Write<Vector> w(v);
                if(r.inside(0)) { v.set(0, 3.0); }
                if(r.inside(1)) { v.set(1, 4.0); }
            }

            double n = norm2(v);
            utopia_test_assert(approxeq(5.0, n));

            n = norm_infty(v);
            utopia_test_assert(approxeq(4.0, n));

            v *= 2.5;
            n = norm2(v);
            utopia_test_assert(approxeq(12.5, n));

            n = norm_infty(v);
            utopia_test_assert(approxeq(10.0, n));
        }

        void quadratic_form()
        {
            const int n = mpi_world_size() * 2;
            Vector x = values(n, 1.);
            Vector b = values(n, 2.);
            Matrix A = values(n, n, 1.);

            double value = 0.5 * dot(x, A * x) + dot(x, b);
            double expected = sum(A) * 0.5 + sum(b);    

            utopia_test_assert(approxeq(value, expected));
            utopia_test_assert(approxeq(value, n*n*0.5 + n*2.));
        }

        void dot_test()
        {
            Vector v1 = zeros(2), v2 = zeros(2);
            {
                auto r = range(v1);
                Write<Vector> w(v1);
                if(r.inside(0)) v1.set(0, 0.0);
                if(r.inside(1)) v1.set(1, 1.0);
            }
            {
                auto r = range(v2);
                Write<Vector> w(v2);
                if(r.inside(0)) v2.set(0, 1.0);
                if(r.inside(1)) v2.set(1, 0.0);
            }

            double v = dot(v1, v2 * 0.1);
            utopia_test_assert(approxeq(0.0, v));
        }

        void dot_product_composition_test()
        {
            Vector v = zeros(2);
            {
                auto r = range(v);
                Write<Vector> w(v);
                if(r.inside(0)) v.set(0, 1.0);
                if(r.inside(1)) v.set(1, 10.0);
            }

            double one = norm2(v) * norm2(v) / dot(v, v);
            utopia_test_assert(approxeq(1.0, one));

            one = norm2(v * (1.0 / Scalar(norm2(v))));
            utopia_test_assert(approxeq(1.0, one));
        }

        void multiply_test()
        {
            if(!is_sparse<Matrix>::value && Traits<Vector>::Backend == PETSC && mpi_world_size() > 1){
                std::cerr << "[Warning] petsc does not support parallel dense matrix-matrix multiplication" << std::endl;
                return;
            }

            Matrix m1 = identity(3, 3);
            {
                Write<Matrix> w(m1, GLOBAL_INSERT);
                m1.set(0, 1, 1);
            }

            Matrix m2 = values(3, 3, 2);
            Matrix m3 = m2 * transpose(m2);
            m3 = transpose(m1) * m3;
            m3 = m2 * m3;
            m3 = m1 * m3;

            each_read(m3, [](SizeType x, SizeType y, double entry) {
                if (x == 0) {
                    utopia_test_assert(entry == 192);
                } else {
                    utopia_test_assert(entry == 96);
                }
            });
        }

        void determinant_test()
        {   
            if(mpi_world_size() > 1) {
                std::cerr << "[Warning] determinant only implemented for serial and small matrices" << std::endl;
                return;
            }

            int n = 3;
            Matrix m = 0.5 * identity(n, n);
            auto expr = det(m);

            double val = expr;
            utopia_test_assert(approxeq(0.125, val));


            Matrix m4 = 0.5 * identity(4, 4);
            double det4 = det(m4);

            utopia_test_assert(approxeq(0.0625, det4));
        }

        void size_test()
        {
            Matrix m = zeros(2, 3);
            Size size = m.size();
            utopia_test_assert(size.get(0) == 2);
            utopia_test_assert(size.get(1) == 3);

            Vector v = zeros(3);
            size = v.size();
            utopia_test_assert(size.get(0) == 3);

            v = m * v;
            size = v.size();
            utopia_test_assert(size.get(0) == 2);
        }

        void binary_min_max()
        {
            const int n = mpi_world_size() * 2;
            Vector one = values(n, 1.);
            Vector two = values(n, 2.);

            Vector actual_min = utopia::min(one, two);
            Vector actual_max = utopia::max(one, two);
            
            utopia_test_assert(approxeq(one, actual_min));
            utopia_test_assert(approxeq(two, actual_max));
        
            actual_min = utopia::min(two, values(n, 1.));
            actual_max = utopia::max(values(n, 2.), one);
        
            utopia_test_assert(approxeq(one, actual_min));
            utopia_test_assert(approxeq(two, actual_max));
        }


        void local_values_test()
        {
            auto k = 15;
            auto m = 4;

            Matrix A = local_values(k, m, 1.);
            Vector x = local_values(k, 1.0);
            Vector x_result = transpose(A)*x; 
            Scalar x_norm = norm2(x_result); 

            utopia_test_assert(x_norm!=0.0);
        }

        void is_subtree()
        {
            Vector v;
            Matrix mat;

            auto expr = -diag(0.1 * v +  mat * v);
            static_assert( !(IsSubTree<Binary<Vector, Vector, Plus>, decltype(expr)>::value), "should be false" );
            static_assert( (IsSubTree<Multiply<Matrix, Vector>, decltype(expr)>::value), "should be true" );
            static_assert( (IsSubTree<Multiply<Matrix, Vector>, decltype(expr)>::value),  "should be true"  );
            static_assert( (IsSubTree<Matrix, decltype(expr)>::value),  "should be true"  );
            static_assert( (IsSubTree<Vector, decltype(expr)>::value),  "should be true"  ); 
        }

        static void print_backend_info()
        {
            if(Utopia::instance().verbose() && mpi_world_rank() == 0) {
                std::cout << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

    public:
        void run()
        {
            print_backend_info();
            UTOPIA_RUN_TEST(is_subtree);
            UTOPIA_RUN_TEST(norm_test);
            UTOPIA_RUN_TEST(dot_test);
            UTOPIA_RUN_TEST(dot_product_composition_test);
            UTOPIA_RUN_TEST(multiply_test);
            UTOPIA_RUN_TEST(determinant_test);
            UTOPIA_RUN_TEST(size_test);
            UTOPIA_RUN_TEST(binary_min_max);
            UTOPIA_RUN_TEST(quadratic_form);
            UTOPIA_RUN_TEST(local_values_test);
        }
    };

    void runAlgebraTest()
    {
        UTOPIA_UNIT_TEST_BEGIN("AlgebraTest");

// #ifdef WITH_BLAS
//         AlgebraTest<Matrixd, Vectord>().run();
// #endif //WITH_BLAS

#ifdef WITH_PETSC
        AlgebraTest<DMatrixd, DVectord>().run();
#endif //WITH_PETSC

// #ifdef WITH_TRILINOS
//         AlgebraTest<TMatrixd, TVectord>().run();
// #endif //WITH_TRILINOS

        UTOPIA_UNIT_TEST_END("AlgebraTest");
    }
}

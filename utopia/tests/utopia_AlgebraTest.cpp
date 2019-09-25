#include "utopia.hpp"
#include "utopia_AlgebraTest.hpp"
#include "utopia_IsSubTree.hpp"

#ifdef WITH_TRILINOS
#include "utopia_trilinos.hpp"
#endif

namespace utopia {

    template<class Matrix, class Vector>
    class SerialAlgebraTest {
    public:
        using Traits = utopia::Traits<Vector>;

        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Expr   = utopia::Multiply<Number<Scalar>, Tensor<Matrix, 2>>;
        using ChosenType = typename ChooseType<Number<Scalar>, Tensor<Matrix, 2>, Number<Scalar>>::Type;
        static const int FT = utopia::Traits<Tensor<Matrix, 2>>::FILL_TYPE;
        static const int Order = Expr::Order;
        using Result = typename TypeAndFill<Traits, Expr>::Type;

        static_assert(Order == 2, "must be a 2nd order tensor");
        static_assert(FT == FillType::DENSE, "must be dense");
        static_assert(std::is_same<ChosenType, Matrix>::value,  "expression must result in matrix");
        static_assert(std::is_same<Result, Matrix>::value, "expression must result in matrix");

        void complicated_test()
        {
            const Scalar lambda = 1.0, mu = 1.0;
            Matrix F = identity(3, 3);
            Matrix H = identity(3, 3);

            Matrix F_inv_t = transpose(inv(F));
            const Scalar J = det(F);
            const Scalar alpha = (1.0 * lambda * std::log(J) - 1.0 * mu);

            const Scalar temp = lambda * inner(F_inv_t, H);
            const Scalar temp2 = inner(F_inv_t, H);
            Matrix mat = mu * H - alpha * F_inv_t * transpose(H) * F_inv_t + lambda * inner(F_inv_t, H) * F_inv_t;
            // Matrix mat = inner(F_inv_t, H) * F_inv_t;

            // std::cout << tree_format((inner(F_inv_t, H) * F_inv_t).get_class()) << std::endl;
        }

        void run()
        {
            UTOPIA_RUN_TEST(complicated_test);
        }

    };

    template<class Vector>
    class VectorAlgebraTest {
    public:
        typedef typename utopia::Traits<Vector>::Scalar Scalar;
        typedef typename utopia::Traits<Vector>::SizeType SizeType;

        static_assert(Traits<Vector>::Order == 1, "Tensor order of vector must be one");
        
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

        void axpy_test()
        {
            int n = 10;
            const Scalar beta = 1.0, omega = 2.0, alpha = 7.0;
            Vector r = values(n, 1.0), 
                   p = values(n, 2.0),
                   v = values(n, 3.0),
                   x = values(n, 500.1),
                   y = values(n, 19.3);

            p = r + beta * (p - omega * v);

            Vector expected = values(n, -3.0);

            utopia_test_assert(approxeq(expected, p));

            Vector s, h;
            s = r - alpha * v;

            expected.set(-20.);
            utopia_test_assert(approxeq(expected, s));

            h = x + alpha * y;
            expected.set(635.2);
            utopia_test_assert(approxeq(expected, h));
        }

        void divide_dots_test()
        {
            int n = 9;
            Vector t = values(n, 1.0), 
                   s = values(n, 2.0);

            const Scalar res = dot(t, s)/dot(t, t);
            utopia_test_assert(approxeq(2.0, res));
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
            UTOPIA_RUN_TEST(norm_test);
            UTOPIA_RUN_TEST(dot_test);
            UTOPIA_RUN_TEST(dot_product_composition_test);
            UTOPIA_RUN_TEST(binary_min_max);
            UTOPIA_RUN_TEST(axpy_test);
            UTOPIA_RUN_TEST(divide_dots_test);
        }

    };

    template<class Matrix, class Vector>
    class AlgebraTest {
    private:
        typedef typename utopia::Traits<Vector>::Scalar Scalar;
        typedef typename utopia::Traits<Vector>::SizeType SizeType;

        static_assert(Traits<Matrix>::Order == 2, "Tensor order of matrix must be 2");

        void nnz_test()
        {
            long n = 100;
            Matrix I = identity(n, n);
            long nnz_I = utopia::nnz(I, 0.);
            utopia_test_assert(nnz_I == n);
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

        void multiply_test()
        {
            if(!is_sparse<Matrix>::value && Traits<Vector>::Backend == PETSC && mpi_world_size() > 1){
                std::cerr << "[Warning] petsc does not support parallel dense matrix-matrix multiplication" << std::endl;
                return;
            }

            Matrix m1 = dense_identity(3, 3);
            {
                Write<Matrix> w(m1, GLOBAL_INSERT);
                m1.c_set(0, 1, 1);
            }

            Matrix m2 = values(3, 3, 2);
            
            Matrix m3 = m2 * transpose(m2);
            //direct variant (1): Matrix m3; m2.transpose_multiply(m2, m3);
            //direct variant (2): Matrix m3; m3.multiply_transpose(m2, m3);

            m3 = transpose(m1) * m3;
            //direct variant: m1.transpose_multiply(m3, m3);

            m3 = m2 * m3;
            //direct variant: m2.multiply(m3, m3);

            m3 = m1 * m3;
            //direct variant:  m1.multiply(m3, m3);

            each_read(m3, [](SizeType i, SizeType /*y*/, double entry) {
                if (i == 0) {
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
            Matrix m = 0.5 * dense_identity(n, n);
            auto expr = det(m);

            double val = expr;
            utopia_test_assert(approxeq(0.125, val));


            Matrix m4 = 0.5 * dense_identity(4, 4);
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
            size = {v.size()};
            utopia_test_assert(size.get(0) == 3);

            v = m * v;
            size = {v.size()};
            utopia_test_assert(size.get(0) == 2);
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

            auto mv = mat * v;
            static_assert( (IsSubTree<decltype(mv), decltype(mv)>::value), "should be true" );

            auto expr = -diag(0.1 * v +  mat * v);
            static_assert( !(IsSubTree<Binary<Vector, Vector, Plus>, decltype(expr)>::value), "should be false" );
            static_assert( (IsSubTree<decltype(mv), decltype(expr)>::value), "should be true" );
            static_assert( (IsSubTree<decltype(mv), decltype(expr)>::value),  "should be true"  );
            static_assert( (IsSubTree<Matrix, decltype(expr)>::value),  "should be true"  );
            static_assert( (IsSubTree<Vector, decltype(expr)>::value),  "should be true"  );
        }

        void trace_test()
        {
            int n = 3;
            Matrix m = 0.5 * dense_identity(n, n);
            Scalar t = trace(m);
            utopia_test_assert(approxeq(t, 1.5, 1e-16));
        }

        void in_place_test()
        {
            int n = 3;
            Matrix oracle = dense_identity(n, n);
            Matrix m = 0.5 * dense_identity(n, n);
            m *= 2.0;

            utopia_test_assert(approxeq(m, oracle, 1e-16));

            oracle *= 1./4.0;
            m /= 4.0;

            utopia_test_assert(approxeq(m, oracle, 1e-16));
        }

        void outer_product_test()
        {
            int n = 5;
            Vector x = values(10, 1.0); 
            Vector y = values(10, 5.0); 


            {
                Write<Vector> rx(x), ry(y);

                each_write(x, [](const SizeType i) -> double {
                    return i+1.0;  }   );

                each_write(y, [](const SizeType i) -> double {
                    return 1./(i+1.);  }   );                
            }            
            
            
            // disp(x); 
            // disp(y); 

            Matrix M = outer(x, y); 
            // disp(M); 

            utopia_test_assert(approxeq(norm2(M), 24.426636618640689, 1e-8));
            
        }



    public:
        void run()
        {
            VectorAlgebraTest<Vector>().run();
            UTOPIA_RUN_TEST(is_subtree);
            UTOPIA_RUN_TEST(multiply_test);
            UTOPIA_RUN_TEST(determinant_test);
            UTOPIA_RUN_TEST(size_test);
            UTOPIA_RUN_TEST(quadratic_form);
            UTOPIA_RUN_TEST(local_values_test);
            UTOPIA_RUN_TEST(nnz_test);
            UTOPIA_RUN_TEST(trace_test);
            UTOPIA_RUN_TEST(in_place_test);
            UTOPIA_RUN_TEST(outer_product_test);
        }
    };

    void runAlgebraTest()
    {
        UTOPIA_UNIT_TEST_BEGIN("AlgebraTest");

#ifdef WITH_BLAS
        AlgebraTest<BlasMatrixd, BlasVectord>().run();
        SerialAlgebraTest<BlasMatrixd, BlasVectord>().run();
#endif //WITH_BLAS

#ifdef WITH_PETSC
        AlgebraTest<PetscMatrix, PetscVector>().run();
#endif //WITH_PETSC

#ifdef WITH_TRILINOS
        VectorAlgebraTest<TpetraVector>().run();
#endif //WITH_TRILINOS

        UTOPIA_UNIT_TEST_END("AlgebraTest");
    }
}

/*
* @Author: Eric Botter
* @Date:   2016-11-07
*/
#include "utopia.hpp"
#include "utopia_AlgebraTest.hpp"
// #include "utopia_Collection.hpp"
#include "utopia_IsSubTree.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class AlgebraTest {
    private:
        typedef typename utopia::Traits<Vector>::Scalar Scalar;

        //FIXME(eric): original norm_test in main.cpp is still there
        void norm_test()
        {
            Vector v = zeros(2);
            {
                Write<Vector> w(v);
                v.set(0, 3.0);
                v.set(1, 4.0);
            }

            double n = norm2(v);
            assert(approxeq(5.0, n));

            n = norm_infty(v);
            assert(approxeq(4.0, n));

            v *= 2.5;
            n = norm2(v);
            assert(approxeq(12.5, n));

            n = norm_infty(v);
            assert(approxeq(10.0, n));
        }

        void dot_test()
        {
            Vector v1 = zeros(2), v2 = zeros(2);
            {
                Write<Vector> w(v1);
                v1.set(0, 0.0);
                v1.set(1, 1.0);
            }
            {
                Write<Vector> w(v2);
                v2.set(0, 1.0);
                v2.set(1, 0.0);
            }

            double v = dot(v1, v2 * 0.1);
            assert(approxeq(0.0, v));
        }

        void dot_product_composition_test()
        {


            Vector v = zeros(2);
            {
                Write<Vector> w(v);
                v.set(0, 1.0);
                v.set(1, 10.0);
            }

            double one = norm2(v) * norm2(v) / dot(v, v);
            assert(approxeq(1.0, one));

            one = norm2(v * (1.0 / Scalar(norm2(v))));
            assert(approxeq(1.0, one));
        }

        void multiply_test()
        {
            if(!is_sparse<Matrix>::value && Traits<Vector>::Backend == PETSC && mpi_world_size() > 1){
                std::cerr << "[Warning] petsc does not support parallel dense matrix-matrix multiplication" << std::endl;
                return;
            }

            Matrix m1 = identity(3, 3);
            {
                Write<Matrix> w(m1);
                m1.set(0, 1, 1);
            }
            Matrix m2 = values(3, 3, 2);
            Matrix m3 = m2 * transpose(m2);
            m3 = transpose(m1) * m3;
            m3 = m2 * m3;
            m3 = m1 * m3;

            each_read(m3, [](SizeType x, SizeType y, double entry) {
                if (x == 0)
                    assert(entry == 192);
                else
                    assert(entry == 96);
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
            assert(approxeq(0.125, val));
        }

        void size_test()
        {
            Matrix m = zeros(2, 3);
            Size size = m.size();
            assert(size.get(0) == 2);
            assert(size.get(1) == 3);

            Vector v = zeros(3);
            size = v.size();
            assert(size.get(0) == 3);

            v = m * v;
            size = v.size();
            assert(size.get(0) == 2);
        }


        void binary_min_max()
        {
            const int n = mpi_world_size() * 2;
            Vector one = values(n, 1.);
            Vector two = values(n, 2.);

            Vector actual_min = utopia::min(one, two);
            Vector actual_max = utopia::max(one, two);
            
            assert(approxeq(one, actual_min));
            assert(approxeq(two, actual_max));
        
            actual_min = utopia::min(two, values(n, 1.));
            actual_max = utopia::max(values(n, 2.), one);
        
            assert(approxeq(one, actual_min));
            assert(approxeq(two, actual_max));
        }

        void vector_selection_test()
        {
            typedef typename utopia::Traits<Vector>::SizeType SizeType;

            const int n = mpi_world_size() * 3;
            Vector v = zeros(n);
            auto r = range(v);

            {
                Write<Vector> w_v(v);
                for(auto i = r.begin(); i < r.end(); ++i) {
                    v.set(i, i);
                }
            }

            std::vector<SizeType> s;
            s.push_back(r.begin());
            s.push_back(r.end() % n);

            Vector selection = v.select(s);
            auto s_r = range(selection);

            {
                Read<Vector> r_s(selection);
                assert(selection.get(s_r.begin()) == r.begin());
                assert(selection.get(s_r.begin() + 1) == (r.end() % n));
            }

            Scalar sum_v_s = sum(v.select(s));
        }

        void matrix_selection_test()
        {
            typedef typename utopia::Traits<Vector>::SizeType SizeType;
            
            const int n = mpi_world_size() * 3;
            Matrix m = zeros(n, n);
            auto rr = row_range(m);

            {
                Write<Matrix> w_m(m);
                for(auto i = rr.begin(); i < rr.end(); ++i) {
                    for(auto j = 0; j < n; ++j) {
                        m.set(i, j, i * n + j);
                    }
                }
            }

            std::vector<SizeType> r_s;
            std::vector<SizeType> c_s;

            r_s.push_back(rr.begin());
            r_s.push_back(rr.begin() + 1);

            c_s.push_back(0);
            c_s.push_back(2);

            Matrix selection = m.select(r_s, c_s);

            {
                auto s_r = row_range(selection);
                Read<Matrix> r_s(selection);
                assert(selection.get(s_r.begin(), 0) == rr.begin() * n);
                assert(selection.get(s_r.begin(), 1) == rr.begin() * n + 2);
            }

            Matrix row_selection = m.select(r_s);

            {
                auto s_r = row_range(row_selection);
                Read<Matrix> r_s(row_selection);

                for(SizeType i = 0; i < n; ++i) {
                    assert(row_selection.get(s_r.begin(), i) == (rr.begin() * n + i));
                }
            }

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


        // void multiply_collections()
        // {
        //     const int n = mpi_world_size() * 3;
        //     Matrix m = 2 * identity(n, n);
        //     auto rr = row_range(m);

        //     Matrix m2 = values(n, n, 0.1);

        //     std::vector<Matrix> matrices(2, m);
        //     std::vector<Matrix> result = m2 * wrap(matrices);
        // }

        static void print_backend_info()
        {
            if(Utopia::Instance().verbose()) {
                std::cout << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

    public:
        void run()
        {
            print_backend_info();
            UTOPIA_RUN_TEST(is_subtree);
            UTOPIA_RUN_TEST(vector_selection_test);
            UTOPIA_RUN_TEST(matrix_selection_test);
            UTOPIA_RUN_TEST(norm_test);
            UTOPIA_RUN_TEST(dot_test);
            UTOPIA_RUN_TEST(dot_product_composition_test);
            UTOPIA_RUN_TEST(multiply_test);
            UTOPIA_RUN_TEST(determinant_test);
            UTOPIA_RUN_TEST(size_test);
            UTOPIA_RUN_TEST(binary_min_max);
            // UTOPIA_RUN_TEST(multiply_collections);
        }
    };

    void runAlgebraTest()
    {
        UTOPIA_UNIT_TEST_BEGIN("AlgebraTest");

        #ifdef WITH_BLAS
            AlgebraTest<Matrixd, Vectord>().run();
        #endif //WITH_BLAS

        #ifdef WITH_PETSC
            AlgebraTest<DMatrixd, DVectord>().run();
        #endif //WITH_PETSC

        UTOPIA_UNIT_TEST_END("AlgebraTest");
    }
}

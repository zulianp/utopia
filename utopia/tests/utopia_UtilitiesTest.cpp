/*
* @Author: Eric Botter
* @Date:   2016-11-11
*/
#include "utopia.hpp"
#include "utopia_AutoDiff.hpp" //simplify_test
#include "utopia_UtilitiesTest.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class UtilitiesTest {
    private:
        void factory_test() {
            Matrix m = identity(2, 2);
            auto size = m.size();
            assert(size.get(0) == 2);
            assert(size.get(1) == 2);

            each_read(m, [](SizeType x, SizeType y, double entry) {
                if (x == y) {
                    assert(approxeq(1.0, entry));
                } else {
                    assert(approxeq(0.0, entry));
                }
            });

            Matrix m2 = values(2, 2, -4);
            size = m2.size();
            assert(size.get(0) == 2);
            assert(size.get(1) == 2);

            each_read(m2, [](SizeType x, SizeType y, double entry) {
                assert(approxeq(-4.0, entry));
            });
        }

        void wrapper_test() {
            Vector v = values(2, 1.0);
            {
                Write<Vector> w(v);
                v.set(1, 2.0);
            }
            Matrix m = values(2, 2, 0.0);
            {
                Write<Matrix> w(m);
                m.set(0, 0, 1.0);
                m.set(1, 1, 1.0);
            }

            Vector res;

            res = m * v * 0.1;

            assert(approxeq(0.1, res.get(0)));
            assert(approxeq(0.2, res.get(1)));
        }

        void range_test() {
            Matrix m1 = identity(3, 3);
            View<Matrix> view = m1.range(0, 1, 0, 3);
            Matrix m2 = view;
            each_read(m2, [](SizeType x, SizeType y, double entry) {
                assert(approxeq(y == 0 ? 1.0 : 0.0, entry));
            });

            #ifdef WITH_PETSC
                //NOTE(eric): range assignment is NYI in PETSc backend
                if (std::is_same<Matrix, DMatrixd>::value) return;
            #endif

            Matrix m3 = m1;
            m3.range(0, 1, 0, 3) = m2;
            each_read(m3, [](SizeType x, SizeType y, double entry) {
                assert(approxeq(x == y ? 1.0 : 0.0, entry));
            });

            Matrix diff = m1 - m3;
            each_read(diff, [](SizeType x, SizeType y, double entry) {
                assert(approxeq(0.0, entry));
            });

            Matrix m4 = values(4, 4, 0.0);
            m4.range(0, 2, 0, 2) = identity(2, 2);
            each_read(m4, [](SizeType x, SizeType y, double entry) {
                assert(approxeq(x == y && x < 2 ? 1.0 : 0.0, entry));
            });
        }

        void factory_and_operations_test()
        {
            const int n = 3;
            Matrix m = zeros(n, n);
            {
                Write<Matrix> w(m);

                m.set(0, 0, 2.0);
                m.set(0, 1, -1.0);

                m.set(1, 0, -1.0);
                m.set(1, 1, 2.0);
                m.set(1, 2, -1.0);

                m.set(2, 1, -1.0);
                m.set(2, 2, -1.0);
            }

            Vector c = (m + 0.1 * identity(n, n)) * values(n, 0.5);
            assert(c.size().get(0) == 3);
        }

        void simplify_test()
        {
            const int n = 2;
            Matrix m = 2.*identity(n, n);
            Vector v = values(n, 1.0);

            // const double ab = 1.0;
            // std::cout << ((ab * m * m + ab * m).getClass()) << std::endl;
            // std::cout << tree_format((ab * m * m + ab * m).getClass()) << std::endl;

            //Useful when applying automatic diff to remove unwanted expressions such as Id and 0
            //For now only works for trees with with certain trees:  Id * (m + 0) * v + 0 *v -> m * v
            //Bug: Id is removed even if it is not in R^(n x n)
            auto expr   = identity(n, n) * (m + zeros(n, n)) * v + zeros(n, n) * v;
            auto s_expr = simplify(expr);

            // disp(tree_format(s_expr.getClass()));

            Vector expected = m * v;
            Vector actual   = s_expr;

            assert(approxeq(expected, actual));
        }

        void variable_test()
        {
            using namespace std;

            auto v = make_shared<Vector>();
            *v = values(20, 1.0);

            auto w = make_shared<Vector>();
            *w = values(20, 3.0);

            Matrix m = identity(20, 20);

            Variable<Vector, 0> var_0(v);
            Variable<Vector, 1> var_1(w);

            auto expr = 0.5 * dot(m * var_0, var_1 - var_0);

            auto &var_0_found = find_variable<Vector, 0>(expr);
            assert(v.get() == &var_0_found.expr());

            auto &var_1_found = find_variable<Vector, 1>(expr);
            assert(w.get() == &var_1_found.expr());

            auto q = make_shared<Vector>();
            *q = values(20, 4.0);

            var_0_found.set(w);
            var_1_found.set(q);

            double result = expr;
            assert(approxeq(10., result));
        }

    public:
        void inline_eval_test()
        {
            int n = 10;
            Vector v = values(n, 1.0);
            Vector res;
            inline_eval(0.1 * v + abs(sqrt(v) - v), res);
            Vector v_exp = values(n, 0.1);
            assert(approxeq(v_exp, res));

            Matrix m     = identity(n, n);
            Matrix m_res = dense(n, n);

            //we will be reading from n to the end of the function
            Read<Matrix> r(m);

            inline_eval(0.1 * m + abs(sqrt(m) - m), m_res);
            each_read(m_res, [](SizeType x, SizeType y, double entry) {
                assert(approxeq(x == y ? 0.1 : 0.0, entry));
            });

            inline_eval((m + m) * v + v, res);
            v_exp = values(n, 3.0);
            assert(approxeq(v_exp, res));

            Number<double> num = 0;
            inline_eval(dot(m*v, 3*v), num);
            assert(approxeq(30.0, num));

            inline_eval(3. * dot(v, v) +
                dot(m * transpose(m),
                    0.1 * m + 0.5 * -0.6*identity(n, n) + values(n, n, 0.0001)
                ), num
            );
            assert(approxeq(28.001, num));
        }

        void run() {
            factory_test();
            wrapper_test();
            range_test();
            factory_and_operations_test();
            simplify_test();
            variable_test();
            inline_eval_test();
        }
    };

    // void wrapper_test_stdvector() {
    //     const std::vector <double> v1{1.0, 2.0, 3.0};
    //     const std::vector <double> v2{1.0, 2.0, 3.0};
    //     std::vector <double> result;
    //
    //     //Wrap and compute
    //     vref(result) = vref(v1) * 0.1 + vref(v2);
    //
    //     auto v_ref_p = vref(v1);
    //     vref(result) = v_ref_p * 0.1 + vref(v2);
    //
    //     const std::vector <double> expected{1.1, 2.2, 3.3};
    //     for (size_t i = 0; i < 3; i++) {
    //         assert(approxeq(result[i], expected[i]));
    //     }
    // }

    void runUtilitiesTest() {
        std::cout << "Begin: UtilitiesTest" << std::endl;

#ifdef WITH_PETSC
        if(mpi_world_size() == 1) {
            UtilitiesTest<DMatrixd, DVectord>().run();
#ifdef WITH_BLAS
            // interoperability
            UtilitiesTest<DMatrixd, Vectord>().inline_eval_test();
            UtilitiesTest<Matrixd, DVectord>().inline_eval_test();
#endif //WITH_BLAS

        } else {
            std::cerr << "[Warning] UtilitiesTest not run for petsc" << std::endl;
        }
#endif //WITH_PETSC

#ifdef WITH_BLAS
            UtilitiesTest<Matrixd, Vectord>().run();
#endif //WITH_BLAS
        // wrapper_test_stdvector();  // doesnt compile on cluster - TODO: fix it

        std::cout << "End:   UtilitiesTest" << std::endl;
    }
}

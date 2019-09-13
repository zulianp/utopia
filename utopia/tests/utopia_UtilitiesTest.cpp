#include "utopia.hpp"
#include "utopia_AutoDiff.hpp" //simplify_test
#include "utopia_UtilitiesTest.hpp"
#include "utopia_Blocks.hpp"
#include "utopia_Eval_Blocks.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class BlockTest {
    public:

        void run() {
            UTOPIA_RUN_TEST(block_test);
        }

    private:
        void block_test()
        {
            int n1 = 5;
            int n2 = 5;
            int n3 = 3;
            int n4 = 4;

            auto id_ptr_11 = std::make_shared<Matrix>(identity(n1, n2));
            auto id_ptr_12 = std::make_shared<Matrix>(2. * identity(n3, n2));
            auto id_ptr_22 = std::make_shared<Matrix>(3. * identity(n3, n4));

            Blocks<Matrix> b_mat(2, 2, {
                id_ptr_11, nullptr,
                id_ptr_12, id_ptr_22
            });

           auto s = size(b_mat);
           utopia_test_assert(s.get(0) == (n1 + n3));
           utopia_test_assert(s.get(1) == (n2 + n4));

           Matrix mat = b_mat;
           Vector ones_1 = values(n2, 2.);
           Vector ones_2 = values(n4, 1.);

           Vector vec = blocks(ones_1, ones_2);
           Vector r = mat * vec;
           Vector r1 = zeros(n1), r2 = zeros(n3);

           undo_blocks(r, r1, r2);

           utopia_test_assert(
            approxeq(
                double(sum(r)),
                2. * size(r1).get(0) + 7. * size(r2).get(0))
            );

           utopia_test_assert(
            approxeq(
                double(sum(r1)),
                2. * size(r1).get(0))
            );

           utopia_test_assert(
            approxeq(
                double(sum(r2)),
                7. * size(r2).get(0))
            );
        }
    };

    template<class Matrix, class Vector>
    class UtilitiesTest {
    private:
        void csv_read_write()
        {
            Path path = Path(Utopia::instance().get("data_path")) / "csv/test.csv";
            CSV csv;

            utopia_test_assert( csv.read(path) );
            utopia_test_assert( csv.write("./out.csv") );

            int val = -1.;
            csv.get(1, 0, val); utopia_test_assert( val == 0 );
            csv.get(1, 1, val); utopia_test_assert( val == 1 );
            csv.get(1, 2, val); utopia_test_assert( val == 2 );

            csv.get(2, 0, val); utopia_test_assert( val == 1 );
            csv.get(2, 1, val); utopia_test_assert( val == 2 );
            csv.get(2, 2, val); utopia_test_assert( val == 3 );
        }

        void factory_test() {
            Matrix m = identity(2, 2);
            auto size = m.size();
            utopia_test_assert(size.get(0) == 2);
            utopia_test_assert(size.get(1) == 2);

            each_read(m, [](SizeType x, SizeType y, double entry) {
                if (x == y) {
                    utopia_test_assert(approxeq(1.0, entry));
                } else {
                    utopia_test_assert(approxeq(0.0, entry));
                }
            });

            Matrix m2 = values(2, 2, -4);
            size = m2.size();
            utopia_test_assert(size.get(0) == 2);
            utopia_test_assert(size.get(1) == 2);

            each_read(m2, [](SizeType, SizeType, double entry) {
                utopia_test_assert(approxeq(-4.0, entry));
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

            {
                Read<Vector> w_res(res);
                utopia_test_assert(approxeq(0.1, res.get(0)));
                utopia_test_assert(approxeq(0.2, res.get(1)));
            }
        }


        //FIXME
        // void range_test() {
        //     Matrix m1 = identity(3, 3);
        //     View<Matrix> view = m1.range(0, 1, 0, 3);
        //     Matrix m2 = view;
        //     each_read(m2, [](SizeType /*x*/, SizeType y, double entry) {
        //         utopia_test_assert(approxeq(y == 0 ? 1.0 : 0.0, entry));
        //     });

        //     #ifdef WITH_PETSC
        //         //NOTE(eric): range assignment is NYI in Petsc backend
        //     if (std::is_same<Matrix, DMatrixd>::value) return;
        //     #endif

        //     Matrix m3 = m1;
        //     m3.range(0, 1, 0, 3) = m2;
        //     each_read(m3, [](SizeType x, SizeType y, double entry) {
        //         utopia_test_assert(approxeq(x == y ? 1.0 : 0.0, entry));
        //     });

        //     Matrix diff = m1 - m3;
        //     each_read(diff, [](SizeType /*x*/, SizeType /*y*/, double entry) {
        //         utopia_test_assert(approxeq(0.0, entry));
        //     });

        //     Matrix m4 = values(4, 4, 0.0);
        //     m4.range(0, 2, 0, 2) = identity(2, 2);
        //     each_read(m4, [](SizeType x, SizeType y, double entry) {
        //         utopia_test_assert(approxeq(x == y && x < 2 ? 1.0 : 0.0, entry));
        //     });
        // }

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
            utopia_test_assert(c.size() == 3);
        }

        //TODO(eric): move this to AutoDiffTest?
        void simplify_test()
        {
            const int n = 2;
            Matrix m = 2.*identity(n, n);
            Vector v = values(n, 1.0);

            // const double ab = 1.0;
            // std::cout << ((ab * m * m + ab * m).getClass()) << std::endl;
            // std::cout << tree_format((ab * m * m + ab * m).getClass()) << std::endl;

            //Useful when applying automatic diff to remove unwanted expressions such as Id and 0
            //For now only works for trees with with certain sub-trees:  Id * (m + 0) * v + 0 *v -> m * v
            //Bug: Id is removed even if it is not in R^(n x n)
            auto expr   = identity(n, n) * (m + zeros(n, n)) * v + zeros(n, n) * v;
            auto s_expr = simplify(expr);

            // disp(tree_format(s_expr.getClass()));

            Vector expected = m * v;
            Vector actual   = s_expr;

            utopia_test_assert(approxeq(expected, actual));
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
            utopia_test_assert(v.get() == &var_0_found.expr());

            auto &var_1_found = find_variable<Vector, 1>(expr);
            utopia_test_assert(w.get() == &var_1_found.expr());

            auto q = make_shared<Vector>();
            *q = values(20, 4.0);

            var_0_found.set(w);
            var_1_found.set(q);

            double result = expr;
            utopia_test_assert(approxeq(10., result));
        }

    public:
        void inline_eval_test()
        {
            int n = 10;
            Vector v = values(n, 1.0);
            Vector res = zeros(n);

            {
                Read<Vector> r_v(v);
                inline_eval(0.1 * v + abs(sqrt(v) - v), res);
            }

            Vector v_exp = values(n, 0.1);
            utopia_test_assert(approxeq(v_exp, res));

            Matrix m     = identity(n, n);
            Matrix m_res = dense(n, n);

            //we will be reading from n to the end of the function
            {
                Read<Matrix> r_m(m);
                inline_eval(0.1 * m + abs(sqrt(m) - m), m_res);
            }

            each_read(m_res, [](SizeType x, SizeType y, double entry) {
                utopia_test_assert(approxeq(x == y ? 0.1 : 0.0, entry));
            });

            {
                Read<Matrix> r_m(m);
                Read<Vector> r_v(v);
                inline_eval((m + m) * v + v, res);
            }

            v_exp = values(n, 3.0);
            utopia_test_assert(approxeq(v_exp, res));

            Number<double> num = 0;

            {
                Read<Matrix> r_m(m);
                Read<Vector> r_v(v);

                inline_eval(dot(m*v, 3*v), num);
            }


            utopia_test_assert(approxeq(30.0, num));

            {
                Read<Matrix> r_m(m);
                Read<Vector> r_v(v);

                inline_eval(3. * dot(v, v) +
                    dot(m * transpose(m),
                        0.1 * m + 0.5 * -0.6*identity(n, n) + values(n, n, 0.0001)
                        ), num
                    );
            }
            utopia_test_assert(approxeq(28.001, num));
        }

        static void print_backend_info()
        {
            mpi_world_barrier();
            if(Utopia::instance().verbose() && mpi_world_rank() == 0) {
                std::cout << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
            mpi_world_barrier();
        }

        void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(csv_read_write);
            UTOPIA_RUN_TEST(factory_test);
            UTOPIA_RUN_TEST(wrapper_test);
            // UTOPIA_RUN_TEST(range_test);
            UTOPIA_RUN_TEST(factory_and_operations_test);
            UTOPIA_RUN_TEST(simplify_test);
            UTOPIA_RUN_TEST(variable_test);
            UTOPIA_RUN_TEST(inline_eval_test);
        }
    };

    void runUtilitiesTest() {
        UTOPIA_UNIT_TEST_BEGIN("UtilitiesTest");

#ifdef WITH_BLAS
        UtilitiesTest<Matrixd, Vectord>().run();
#endif //WITH_BLAS

#ifdef WITH_PETSC
        BlockTest<DSMatrixd, DVectord>().run();


        if(mpi_world_size() == 1) {
            // UtilitiesTest<DMatrixd, DVectord>().run(); //FIXME
            BlockTest<DMatrixd, DVectord>().run();
#ifdef WITH_BLAS
            // interoperability
            // UtilitiesTest<DMatrixd, Vectord>().inline_eval_test();
            // UtilitiesTest<Matrixd, DVectord>().inline_eval_test();
#endif //WITH_BLAS

        }

#endif //WITH_PETSC

//         if(mpi_world_size() == 1) {
// #ifdef WITH_TRILINOS
//             BlockTest<TSMatrixd, TVectord>().run();
// #endif //WITH_TRILINOS
//         }

        UTOPIA_UNIT_TEST_END("UtilitiesTest");
    }
}

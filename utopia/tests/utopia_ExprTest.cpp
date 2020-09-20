#include <utility>
#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_Eval_Residual.hpp"
#include "utopia_For.hpp"
#include "utopia_Testing.hpp"

namespace utopia {

    template <typename Matrix, typename Vector>
    class ExpressionTests {
    public:
        using Traits = utopia::Traits<Vector>;
        using SizeType = typename Traits::SizeType;
        using Scalar = typename Traits::Scalar;
        using Layout = typename Traits::Layout;
        using MatrixLayout = typename Traits::MatrixLayout;

        using Comm = typename Traits::Communicator;

        ExpressionTests(const Comm &comm, const SizeType &n) {
            vec_layout = layout(comm, n, n * comm.size());
            mat_layout = square_matrix_layout(vec_layout);
        }

        void run() {
            UTOPIA_RUN_TEST(axpy_test);
            UTOPIA_RUN_TEST(e_div_test);
            UTOPIA_RUN_TEST(e_mul_test);
            UTOPIA_RUN_TEST(for_each_loop_test);
            UTOPIA_RUN_TEST(multi_reduce_test);
            UTOPIA_RUN_TEST(mv_test);
            UTOPIA_RUN_TEST(negate_alpha_test);
            UTOPIA_RUN_TEST(negate_test);
            UTOPIA_RUN_TEST(parallel_each_write_test);
            UTOPIA_RUN_TEST(quad_form_test);
            UTOPIA_RUN_TEST(residual_test);
            UTOPIA_RUN_TEST(transform_test);
            UTOPIA_RUN_TEST(mat_copy);
            UTOPIA_RUN_TEST(reciprocal_test);
            UTOPIA_RUN_TEST(max_min_test);
            UTOPIA_RUN_TEST(multi_axpy);
            UTOPIA_RUN_TEST(inv_diag);
            UTOPIA_RUN_TEST(comp_mat);
            UTOPIA_RUN_TEST(bratu_grad);
            UTOPIA_RUN_TEST(diag_mult);
            UTOPIA_RUN_TEST(rotate_test);
            UTOPIA_RUN_TEST(e_pseudo_inv_test);
            UTOPIA_RUN_TEST(mat_vec_multiply_test);
            UTOPIA_RUN_TEST(vec_add_add_add_test);
            UTOPIA_RUN_TEST(convert_test);
            UTOPIA_RUN_TEST(emul_test);
            UTOPIA_RUN_TEST(mat_transp_mult_test);
        }

        // void convert_test()
        // {
        //     #ifdef UTOPIA_WITH_PETSC
        //         if(Traits<Vector>::Backend == PETSC) {
        //             Vector ut_vec = local_zeros(10);

        //             Vec x;
        //             VecDuplicate(raw_type(ut_vec), &x);
        //             convert(ut_vec, x);
        //             ut_vec *= 500;

        //             UTOPIA_NO_ALLOC_BEGIN("Ellipse::convert_test");
        //             convert(x, ut_vec);
        //             UTOPIA_NO_ALLOC_END();
        //         }
        //     #endif
        // }
        void convert_test() {
            // FIXME create converts for Trilinos then enable this
            // Vector old_vec = local_zeros(10);
            // Vector ut_vec  = local_zeros(10);

            // convert(ut_vec, raw_type(old_vec));
            // ut_vec *= 500;

            // UTOPIA_NO_ALLOC_BEGIN("convert_test");
            // convert(raw_type(old_vec), ut_vec);
            // UTOPIA_NO_ALLOC_END();
        }

        void vec_add_add_add_test() {
            Vector x(vec_layout, 1.0);
            Vector y(vec_layout, 2.0);
            Vector z(vec_layout, 3.0);

            UTOPIA_NO_ALLOC_BEGIN("vec_add_add_add_test");
            x += y + z;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("vec_add_add_add_test2");
            x = x + y + z;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("vec_add_add_add_test3");
            z = x + y + z;
            UTOPIA_NO_ALLOC_END();
        }

        void mat_transp_mult_test() {
            Matrix H;
            H.sparse(mat_layout, 3, 2);
            assemble_laplacian_1D(H);
            Matrix D = diag(diag(H));
            H = H + transpose(H) - D;
        }

        void mat_vec_multiply_test() {
            Matrix H;
            H.sparse(mat_layout, 3, 2);
            assemble_laplacian_1D(H);

            Vector x(vec_layout, 1.0);
            Vector y(vec_layout, 1.0);

            // UTOPIA_NO_ALLOC_BEGIN("mat_vec_multiply_test1");
            // x = H*x;
            // UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("mat_vec_multiply_test2");
            y = H * x;
            UTOPIA_NO_ALLOC_END();
        }

        void emul_test() {
            Vector x(vec_layout, 1.0);
            Vector y(vec_layout, 10.0);
            Vector z(vec_layout, 10.0);

            double scaling_param = 10.0;

            UTOPIA_NO_ALLOC_BEGIN("emul_test");
            z = e_mul(x, y);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("emul2_test");
            y = scaling_param * e_mul(x, y);
            UTOPIA_NO_ALLOC_END();
        }

        void negate_test() {
            Vector x(vec_layout, 1.0);
            Vector y(vec_layout, 10.0);

            UTOPIA_NO_ALLOC_BEGIN("negate_test");
            x = -x;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("copy_test");
            x = y;
            UTOPIA_NO_ALLOC_END();
        }

        void mv_test() {
            Vector x(vec_layout, 1.0);
            Vector b(vec_layout, 1.0);
            Vector p(vec_layout, 1.0);

            Matrix A;
            A.sparse(mat_layout, 3, 3);
            assemble_laplacian_1D(A);

            UTOPIA_NO_ALLOC_BEGIN("mv_test1");
            p = A * x + b;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("mv_test2");
            p = A * x - b;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("mv_test3");
            p = b - A * x;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("mv_test4");
            p = b + A * x;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("mv_test5");
            p = 0.5 * b + A * x;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("mv_test6");
            p = 0.5 * b - A * x;
            UTOPIA_NO_ALLOC_END();

            // This will always create a copy (A * x needs to put the result somewhere, unless the backend
            // specialization is used and does not allocate anything) a better way would be to use a third vector p = b
            // - A * x;
            b = b - A * x;

            UTOPIA_NO_ALLOC_BEGIN("mv_test7");
            p = b - A * x;
            UTOPIA_NO_ALLOC_END();
        }

        void mat_copy() {
            Matrix I;
            I.identity(mat_layout);
            Matrix D;
            D.identity(mat_layout);
            Vector v = diag(I);

            // I do not know how relevant are these tests, as sparsity pattern might be different...
            // create same_sparsity_copy()
            // UTOPIA_NO_ALLOC_BEGIN("mat_copy1");
            // FIME
            // D = I;
            // UTOPIA_NO_ALLOC_END();

            // UTOPIA_NO_ALLOC_BEGIN("mat_copy2");
            // FIXME this it is equivalent to making a copy (but it could be treated as an AXPY)
            D = -1.0 * I;
            D = -I;
            // UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("mat_copy3");
            D -= I;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("mat_copy4");
            D += I;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("mat_copy5");
            D += 5.0 * D;
            UTOPIA_NO_ALLOC_END();

            D += Matrix(diag(v));

            // The next lines seem to crash when using Trilinos/Code coverage info/specific version of Kokkos/Tpetra
            // #ifndef UTOPIA_WITH_CODE_COVERAGE

            UTOPIA_NO_ALLOC_BEGIN("mat_copy5");  //
            // https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatDiagonalSet.html#MatDiagonalSet
            D.shift_diag(v);
            D += diag(v);  // bad way D += Matrix(diag(v));
            UTOPIA_NO_ALLOC_END();

            // UTOPIA_NO_ALLOC_BEGIN("mat_copy6");
            // FIME still creates a temporary (but now it is just a vector)
            D -= diag(v);  // bad way D -= Matrix(diag(v));
                           // UTOPIA_NO_ALLOC_END();

            // #else
            // #warning "FIXME: D += diag(v) test si commented out because it generates errors with code-coverage
            // instrumentation" #endif
        }

        void e_mul_test() {
            Vector x(vec_layout, 1.0);
            Vector y(vec_layout, 2.0);
            Vector z(vec_layout, 0.0);

            Matrix D = diag(z);

            UTOPIA_NO_ALLOC_BEGIN("e_mul_test");
            z = e_mul(x, y);
            UTOPIA_NO_ALLOC_END();
            utopia_test_assert(approxeq(sum(z), 2 * x.size()));
        }

        void max_min_test() {
            Vector x(vec_layout, 1.0);
            Vector y(vec_layout, 2.0);
            Vector z(vec_layout, 0.0);

            UTOPIA_NO_ALLOC_BEGIN("max_min_test1");
            z = utopia::max(x, y);
            UTOPIA_NO_ALLOC_END();
            utopia_test_assert(approxeq(sum(z), 2 * x.size()));

            UTOPIA_NO_ALLOC_BEGIN("max_min_test2");
            z = utopia::min(x, y);
            UTOPIA_NO_ALLOC_END();
            utopia_test_assert(approxeq(sum(z), x.size()));

            UTOPIA_NO_ALLOC_BEGIN("max_min_test3");
            z = utopia::max(utopia::min(x, y), y);
            UTOPIA_NO_ALLOC_END();

            utopia_test_assert(approxeq(sum(z), 2 * x.size()));
        }

        void reciprocal_test() {
            Vector x(vec_layout, 2.0);

            UTOPIA_NO_ALLOC_BEGIN("reciprocal_test");
            x = 1. / x;
            UTOPIA_NO_ALLOC_END();
            utopia_test_assert(approxeq(sum(x), 0.5 * x.size()));
        }

        void e_div_test() {
            Vector x(vec_layout, 6.0);
            Vector z(vec_layout, 3.0);

            UTOPIA_NO_ALLOC_BEGIN("e_div_test");
            z = x / z;

            Scalar sum_z = sum(z);
            utopia_test_assert(approxeq(sum_z, 2.0 * x.size()));
            z = x / x;
            sum_z = sum(z);
            utopia_test_assert(approxeq(sum_z, 1.0 * x.size()));

            z = z / x;
            sum_z = sum(z);
            utopia_test_assert(approxeq(sum_z, 1.0 / 6.0 * x.size()));

            UTOPIA_NO_ALLOC_END();
        }

        void transform_test() {
            Vector x(vec_layout, 1.0);

            auto x_view = view_device(x);
            parallel_for(
                range_device(x), UTOPIA_LAMBDA(const SizeType &i) { x_view.set(i, (i + 1) * x_view.get(i)); });

            Scalar expected = ((x.size() + 1) * x.size()) / 2.0;
            Scalar sum_x = sum(x);

            utopia_test_assert(approxeq(sum_x, expected, 1e-10));
        }

        void negate_alpha_test() {
            Vector x(vec_layout, 1.0);
            Vector y(vec_layout, 1.0);

            UTOPIA_NO_ALLOC_BEGIN("negate_alpha_test");
            y = -0.5 * x;
            UTOPIA_NO_ALLOC_END();
        }

        void quad_form_test() {
            Vector x(vec_layout, 1.0);
            Vector y(vec_layout, 2.0);

            auto expr = 0.5 * dot(x, y) - 0.5 * dot(x, y);

            Scalar val = expr;

            utopia_test_assert(approxeq(val, 0.0));
        }

        void residual_test() {
            Vector x(vec_layout, 1.0);
            Vector b(vec_layout, 1.0);
            Matrix A;
            A.sparse(mat_layout, 3, 3);
            Vector r(vec_layout, 0.0);

            UTOPIA_NO_ALLOC_BEGIN("residual_test");
            r = x - b;
            UTOPIA_NO_ALLOC_END();
        }

        void axpy_test() {
            Vector x(vec_layout, 1.0);
            Vector y(vec_layout, 1.0);
            Vector p(vec_layout, 2.0);

            UTOPIA_NO_ALLOC_BEGIN("axpy_test0");
            x = x - 0.5 * p;
            y = x - 0.5 * p;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("axpy_test1");
            x.set(1.0);
            x = x + 0.5 * p;
            UTOPIA_NO_ALLOC_END();
            Scalar val = sum(x);
            utopia_test_assert(approxeq(val, x.size() * 2.0));

            UTOPIA_NO_ALLOC_BEGIN("axpy_test2");
            x.set(2.0);
            y.set(1.0);
            x = y + 0.5 * x;
            UTOPIA_NO_ALLOC_END();
            val = sum(x);
            utopia_test_assert(approxeq(val, x.size() * 2.0));
        }

        void for_each_loop_test() {
            Vector x(vec_layout, 2.);
            Vector y(vec_layout, 1.);
            Vector z(vec_layout, 0.);

            {
                auto d_x = const_local_view_device(x);
                auto d_y = const_local_view_device(y);
                auto d_z = local_view_device(z);

                parallel_for(
                    local_range_device(z), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar xi = d_x.get(i);
                        const Scalar yi = d_y.get(i);
                        d_z.set(i, xi - yi);
                    });
            }

            Scalar sum_z = sum(z);
            utopia_test_assert(approxeq(Scalar(x.size()), sum_z));
        }

        void parallel_each_write_test() {
            Vector x(vec_layout, 2.);
            Vector y(vec_layout, 1.);
            Vector z(vec_layout, 0.);

            {
                auto d_x = const_local_view_device(x);
                auto d_y = const_local_view_device(y);
                auto d_z = local_view_device(z);

                parallel_for(
                    local_range_device(z), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar xi = d_x.get(i);
                        const Scalar yi = d_y.get(i);
                        d_z.set(i, xi - yi);
                    });
            }

            Scalar sum_z = sum(z);
            utopia_test_assert(approxeq(Scalar(x.size()), sum_z));
        }

        void multi_reduce_test() {
            Vector x(vec_layout, 2.);
            Vector y(vec_layout, 1.);

            {
                auto x_view = view_device(x);
                auto y_view = view_device(y);

                parallel_for(
                    range_device(x), UTOPIA_LAMBDA(const SizeType &i) {
                        x_view.set(i, -(i + 1.0));
                        y_view.set(i, (i + 1.0));
                    });
            }

            const Scalar m = multi_min(x, y);
            utopia_test_assert(approxeq(m, Scalar(-x.size())));
        }

        void multi_axpy() {
            // Assign<Vec, Minus<Plus<Vec, Multiplies<Number, Vec>>, Multiplies<Number, Vec>>>
            Vector a(vec_layout, 2.);
            Vector b(vec_layout, 2.);
            Vector c(vec_layout, 2.);
            Vector result(vec_layout, 0.0);

            const Scalar alpha = 1.0, beta = 2.0;

            UTOPIA_NO_ALLOC_BEGIN("multi_axpy_1");
            result = a + (alpha * b) - (beta * c);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("multi_axpy_2");
            c = a + (alpha * b) - (beta * c);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("multi_axpy_3");
            b = a + (alpha * b) - (beta * c);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("multi_axpy_4");
            b = a - (alpha * b) + (beta * c);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("multi_axpy_5");
            c = a - (alpha * b) + (beta * c);
            UTOPIA_NO_ALLOC_END();
        }

        void inv_diag() {
            Matrix H;
            H.sparse(mat_layout, 3, 2);
            Vector d(vec_layout, 0.0);

            assemble_laplacian_1D(H);

            UTOPIA_NO_ALLOC_BEGIN("inv_diag");
            d = 1. / diag(H);
            UTOPIA_NO_ALLOC_END();
        }

        void comp_mat() {
            Matrix H1;
            H1.identity(mat_layout);
            Matrix H2;
            H2.identity(mat_layout, 2.0);
            Matrix H3;
            H3.identity(mat_layout, 3.0);
            Matrix R = H1;

            if (Traits::Backend == PETSC) {
                UTOPIA_NO_ALLOC_BEGIN("comp_mat");
                R = H3 + H1 * H2;
                UTOPIA_NO_ALLOC_END();
            } else {
                // FIXME
                R = H3 + H1 * H2;
            }

            Matrix Id;
            Id.identity(mat_layout, 5.0);

            utopia_test_asserteq(Id, R, 0.0);
        }

        void bratu_grad() {
            Vector result(vec_layout, 0.0);
            Vector x(vec_layout, 0.0);

            UTOPIA_NO_ALLOC_BEGIN("bratu_grad");
            result = x - (0.5 * exp(x));
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("bratu_grad2");
            result = x - (0.5 * 0.2 * exp(x));
            UTOPIA_NO_ALLOC_END();
        }

        void diag_mult() {
            // FIXME
            if (Traits::Backend == PETSC) {
                Vector d(vec_layout, 2);
                Matrix B;
                B.identity(mat_layout, 2);
                Matrix T;
                T.sparse(mat_layout, 1, 0);

                // UTOPIA_NO_ALLOC_BEGIN("diag_mult");
                T = diag(1. / d) * B;
                // UTOPIA_NO_ALLOC_END();

                Matrix Id;
                Id.identity(mat_layout);
                utopia_test_assert(approxeq(Id, T));
            }
        }

        void rotate_test() {
            Vector x1(vec_layout, 1.0);
            Vector x2(vec_layout, 2.0);
            Vector x3(vec_layout, 3.0);
            Vector x4(vec_layout, 3.0);

            std::vector<Vector> vecs(4);
            vecs[0] = x1;
            vecs[1] = x2;
            vecs[2] = x3;
            vecs[3] = x4;

            // testing that move operations are implemented
            UTOPIA_NO_ALLOC_BEGIN("rotate_test");
            std::rotate(vecs.begin(), vecs.begin() + 1, vecs.end());
            UTOPIA_NO_ALLOC_END();
        }

        void e_pseudo_inv_test() {
            Vector x1, x2;

            x1.values(vec_layout, 1.0);
            e_pseudo_inv(x1, x2);

            utopia_test_assert(!x2.empty());
            utopia_test_assert(range(x1) == range(x2));

            UTOPIA_NO_ALLOC_BEGIN("e_pseudo_inv_test");
            e_pseudo_inv(x1, x2);
            UTOPIA_NO_ALLOC_END();
        }

        Layout vec_layout;
        MatrixLayout mat_layout;
    };

    void expr() {
#ifdef UTOPIA_WITH_PETSC
        ExpressionTests<PetscMatrix, PetscVector>(PetscCommunicator(), 10).run();
#endif

#ifdef UTOPIA_WITH_TRILINOS
        ExpressionTests<TpetraMatrixd, TpetraVectord>(TrilinosCommunicator(), 10).run();
#endif
    }

    UTOPIA_REGISTER_TEST_FUNCTION(expr);
}  // namespace utopia

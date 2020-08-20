#include "utopia.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_Assert.hpp"
#include "utopia_IsSubTree.hpp"
#include "utopia_Testing.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#ifdef WITH_TRILINOS
#include "utopia_trilinos.hpp"
#endif

namespace utopia {

    template <class Matrix, class Vector>
    class SerialAlgebraTest {
    public:
        using Traits = utopia::Traits<Vector>;

        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Expr = utopia::Multiply<Number<Scalar>, Tensor<Matrix, 2>>;
        using ChosenType = typename ChooseType<Number<Scalar>, Tensor<Matrix, 2>, Number<Scalar>>::Type;
        static const int FT = utopia::Traits<Tensor<Matrix, 2>>::FILL_TYPE;
        static const int Order = Expr::Order;
        using Result = typename TypeAndFill<Traits, Expr>::Type;

        static_assert(Order == 2, "must be a 2nd order tensor");
        static_assert(FT == FillType::DENSE, "must be dense");
        static_assert(std::is_same<ChosenType, Matrix>::value, "expression must result in matrix");
        static_assert(std::is_same<Result, Matrix>::value, "expression must result in matrix");

        void complicated_test() {
            const Scalar lambda = 1.0, mu = 1.0;

            auto lo = serial_layout(3, 3);
            Matrix F, H;
            F.identity(lo, 1.0);
            H.identity(lo, 1.0);

            Matrix F_inv_t = transpose(inv(F));
            const Scalar J = det(F);
            const Scalar alpha = (1.0 * lambda * std::log(J) - 1.0 * mu);

            Matrix mat = mu * H - alpha * F_inv_t * transpose(H) * F_inv_t + lambda * inner(F_inv_t, H) * F_inv_t;
            // utopia::out() <<tree_format((inner(F_inv_t, H) * F_inv_t).get_class()) << std::endl;
            utopia_test_assert(SizeType(3) == mat.rows());
            utopia_test_assert(SizeType(3) == mat.cols());
        }

        void run() { UTOPIA_RUN_TEST(complicated_test); }
    };

    template <class Vector>
    class VectorAlgebraTest {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        static_assert(Traits::Order == 1, "Tensor order of vector must be one");

        Comm world;

        void norm_test() {
            Vector v(serial_layout(2), 0.0);

            {
                auto r = range(v);
                Write<Vector> w(v);
                if (r.inside(0)) {
                    v.set(0, 3.0);
                }
                if (r.inside(1)) {
                    v.set(1, 4.0);
                }
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

        void dot_test() {
            auto lo = serial_layout(2);

            Vector v1(lo, 0.0), v2(lo, 0.0);
            {
                auto r = range(v1);
                Write<Vector> w(v1);
                if (r.inside(0)) {
                    v1.set(0, 0.0);
                }
                if (r.inside(1)) {
                    v1.set(1, 1.0);
                }
            }
            {
                auto r = range(v2);
                Write<Vector> w(v2);
                if (r.inside(0)) {
                    v2.set(0, 1.0);
                }
                if (r.inside(1)) {
                    v2.set(1, 0.0);
                }
            }

            double v = dot(v1, v2 * 0.1);
            utopia_test_assert(approxeq(0.0, v));
        }

        void dot_product_composition_test() {
            auto lo = serial_layout(2);

            Vector v(lo, 0.0);
            {
                auto r = range(v);
                Write<Vector> w(v);
                if (r.inside(0)) {
                    v.set(0, 1.0);
                }
                if (r.inside(1)) {
                    v.set(1, 10.0);
                }
            }

            double one = norm2(v) * norm2(v) / dot(v, v);
            utopia_test_assert(approxeq(1.0, one));

            one = norm2(v * (1.0 / Scalar(norm2(v))));
            utopia_test_assert(approxeq(1.0, one));
        }

        void reduce_test() {
            auto lo = serial_layout(5);
            Vector x(lo, 1.0);

            Scalar min_val = min(x);
            Scalar max_val = max(x);
            Scalar sum_val = sum(x);

            utopia_test_assert(approxeq(1.0, min_val));
            utopia_test_assert(approxeq(1.0, max_val));
            utopia_test_assert(approxeq(5.0, sum_val));

            x.set(-1.0);

            min_val = min(x);
            max_val = max(x);
            sum_val = sum(x);

            utopia_test_assert(approxeq(-1.0, min_val));
            utopia_test_assert(approxeq(-1.0, max_val));
            utopia_test_assert(approxeq(-5.0, sum_val));
        }

        void binary_min_max() {
            const int n = mpi_world_size() * 2;
            auto lo = layout(world, 2, n);

            Vector one(lo, 1.);
            Vector two(lo, 2.);

            Vector actual_min = utopia::min(one, two);
            Vector actual_max = utopia::max(one, two);

            utopia_test_assert(approxeq(one, actual_min));
            utopia_test_assert(approxeq(two, actual_max));

            // FIXME
            // actual_min = utopia::e_min(two, 1.);
            // actual_max = utopia::e_max(2., one);

            actual_min = two;
            actual_min.e_min(1.0);

            actual_max = one;
            actual_max.e_max(2.0);

            utopia_test_assert(approxeq(one, actual_min));
            utopia_test_assert(approxeq(two, actual_max));
        }

        void axpy_test() {
            const SizeType n = 10;
            auto v_layout = layout(world, Traits::decide(), n);

            const Scalar beta = 1.0, omega = 2.0, alpha = 7.0;
            Vector r(v_layout, 1.0), p(v_layout, 2.0), v(v_layout, 3.0), x(v_layout, 500.1), y(v_layout, 19.3);

            p = r + beta * (p - omega * v);

            Vector expected(v_layout, -3.0);

            utopia_test_assert(approxeq(expected, p));

            Vector s, h;
            s = r - alpha * v;

            expected.set(-20.);
            utopia_test_assert(approxeq(expected, s));

            h = x + alpha * y;
            expected.set(635.2);
            utopia_test_assert(approxeq(expected, h));

            h = x + (y * alpha);
            utopia_test_assert(approxeq(expected, h));

            h.set(1.0);
            x.set(1.0);

            h = x + (h * alpha);
            expected.set(1.0 + alpha);
            utopia_test_assert(approxeq(expected, h));
        }

        void divide_dots_test() {
            SizeType n = 9;

            auto v_layout = layout(world, Traits::decide(), n);
            Vector t(v_layout, 1.0), s(v_layout, 2.0);

            const Scalar res = dot(t, s) / dot(t, t);
            utopia_test_assert(approxeq(2.0, res));
        }

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

    public:
        void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(norm_test);
            UTOPIA_RUN_TEST(dot_test);
            UTOPIA_RUN_TEST(dot_product_composition_test);
            UTOPIA_RUN_TEST(reduce_test);
            UTOPIA_RUN_TEST(binary_min_max);
            UTOPIA_RUN_TEST(axpy_test);
            UTOPIA_RUN_TEST(divide_dots_test);
        }
    };

    template <class Matrix, class Vector>
    class SparseAlgebraTest {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        static_assert(Traits::Order == 2, "Tensor order of matrix must be 2");

        SparseAlgebraTest() : n(is_dense<Matrix>::value ? 600 : 8000) {}

        Comm world;
        const int n;

        void sparse_chop_test() {
            Matrix M;

            M.sparse(layout(world, Traits::decide(), Traits::decide(), n, n), 3, 3);

            assemble_laplacian_1D(M);

            chop_smaller_than(M, 1e-13);
            chop_greater_than(M, 1e-13);

            utopia_test_assert(approxeq(0.0, norm2(M), 1e-13));
        }

        void transform_test() {
            Matrix M;

            M.sparse(layout(world, Traits::decide(), Traits::decide(), n, n), 3, 3);

            assemble_laplacian_1D(M);
            Matrix M_abs = abs(M);
            const Scalar sum_M_abs = sum(M_abs);

            Scalar expected = 2 * 2 + (n - 2) * 4;
            utopia_test_assert(approxeq(expected, sum_M_abs, 1e-13));

            Matrix sqrt_M_abs = sqrt(M_abs);
            const Scalar sum_sqrt_M_abs = sum(sqrt_M_abs);

            expected = 2 * 2 + (n - 2) * (2 + std::sqrt(2.0));
            utopia_test_assert(approxeq(expected, sum_sqrt_M_abs, 1e-8));

            Matrix M_abs_2 = pow2(sqrt_M_abs);
            utopia_test_assert(approxeq(M_abs, M_abs_2, 1e-8));
        }

        void transform_ijv_test() {
            Matrix M;

            M.sparse(layout(world, Traits::decide(), Traits::decide(), 3, 3), 3, 3);

            assemble_laplacian_1D(M);

            M *= 0.0;

            assemble_laplacian_1D(M);

            auto rr = M.row_range();
            auto nc = M.cols();

            M.transform_ijv(UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &v)->Scalar {
                utopia_test_assert(rr.inside(i));
                utopia_test_assert(j < nc);
                return v;
            });

            M.read(UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &v) {
                utopia_test_assert(rr.inside(i));
                utopia_test_assert(j < nc);
                utopia_test_assert(v >= -1.0);
                utopia_test_assert(v <= 2.0);
            });
        }

        void transpose_test() {
            Matrix M, M_copy;

            M.sparse(layout(world, Traits::decide(), Traits::decide(), 3, 3), 3, 3);

            assemble_laplacian_1D(M);
            M_copy.copy(M);

            M = transpose(M);

            utopia_test_asserteq(M, M_copy, device::epsilon<Scalar>());
        }

        void run() {
            UTOPIA_RUN_TEST(sparse_chop_test);
            UTOPIA_RUN_TEST(transform_test);
            UTOPIA_RUN_TEST(transform_ijv_test);
            UTOPIA_RUN_TEST(transpose_test);
        }
    };

    template <class Matrix, class Vector>
    class DenseAlgebraTest {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        static_assert(Traits::Order == 2, "Tensor order of matrix must be 2");

        Comm world;

        void nnz_test() {
            long n = 100;
            Matrix I;

            I.identity(layout(world, Traits::decide(), Traits::decide(), n, n));

            long nnz_I = utopia::nnz(I, 0.);
            utopia_test_assert(nnz_I == n);
        }

        void quadratic_form() {
            const int n = world.size() * 2;
            auto vec_layout = layout(world, 2, n);
            auto mat_layout = layout(world, 2, 2, n, n);

            Vector x(vec_layout, 1.);
            Vector b(vec_layout, 2.);
            Matrix A;

            A.dense(mat_layout, 1.0);

            double value = 0.5 * dot(x, A * x) + dot(x, b);
            double expected = sum(A) * 0.5 + sum(b);

            utopia_test_assert(approxeq(value, expected));
            utopia_test_assert(approxeq(value, n * n * 0.5 + n * 2.));
        }

        void multiply_test() {
            // if(!is_sparse<Matrix>::value && Traits<Vector>::Backend == PETSC && mpi_world_size() > 1){
            if (Traits::Backend == PETSC && world.size() > 1) {
                std::cerr << "[Warning] petsc does not support parallel dense matrix-matrix multiplication"
                          << std::endl;
                return;
            }

            if (world.size() > 3) {
                std::cerr << "[Warning] multiply_test not run for comm size > 3" << std::endl;
                return;
            }

            const SizeType n = 3;
            auto mat_layout = layout(world, Traits::decide(), Traits::decide(), n, n);

            Matrix m1, m2, m3;
            m1.dense_identity(mat_layout);

            {
                Write<Matrix> w(m1, GLOBAL_INSERT);
                m1.c_set(0, 1, 1);
            }

            m2.dense(mat_layout);
            m2.set(2.0);

            m3 = m2 * transpose(m2);
            // direct variant (1): Matrix m3; m2.transpose_multiply(m2, m3);
            // direct variant (2): Matrix m3; m3.multiply_transpose(m2, m3);

            m3 = transpose(m1) * m3;
            // direct variant: m1.transpose_multiply(m3, m3);

            m3 = m2 * m3;
            // direct variant: m2.multiply(m3, m3);

            m3 = m1 * m3;
            // direct variant:  m1.multiply(m3, m3);

            m3.read([](SizeType i, SizeType /*y*/, double entry) {
                if (i == 0) {
                    utopia_test_assert(entry == 192);
                } else {
                    utopia_test_assert(entry == 96);
                }
            });
        }

        void determinant_test() {
            Matrix m, m4;

            SizeType n = 3;
            m.dense_identity(serial_layout(n, n), 0.5);

            utopia_test_assert(m.rows() == n);
            utopia_test_assert(m.cols() == n);

            auto expr = det(m);

            double val = expr;
            utopia_test_assert(approxeq(0.125, val));

            n = 4;
            m4.dense_identity(serial_layout(n, n), 0.5);

            utopia_test_assert(m4.rows() == n);
            utopia_test_assert(m4.cols() == n);

            double det4 = det(m4);

            utopia_test_assert(approxeq(0.0625, det4));
        }

        void size_test() {
            Matrix m(serial_layout(2, 3));

            Size size = m.size();
            utopia_test_assert(size.get(0) == 2);
            utopia_test_assert(size.get(1) == 3);

            Vector v(serial_layout(3));
            size = {v.size()};
            utopia_test_assert(size.get(0) == 3);

            v = m * v;
            size = {v.size()};
            utopia_test_assert(size.get(0) == 2);
        }

        void local_values_test() {
            const SizeType k = 15;
            const SizeType m = 4;
            const SizeType global_k = k * world.size();
            const SizeType global_m = m * world.size();

            Matrix A;
            Vector x, x_result;

            A.dense(layout(world, k, m, global_k, global_m));
            A.set(1.0);
            x.values(layout(world, k, global_k), 1.0);
            x_result = transpose(A) * x;

            Scalar x_norm = norm2(x_result);
            utopia_test_assert(x_norm != 0.0);
        }

        void is_subtree() {
            Vector v;
            Matrix mat;

            auto mv = mat * v;
            static_assert((IsSubTree<decltype(mv), decltype(mv)>::value), "should be true");

            auto expr = -diag(0.1 * v + mat * v);
            static_assert(!(IsSubTree<Binary<Vector, Vector, Plus>, decltype(expr)>::value), "should be false");
            static_assert((IsSubTree<decltype(mv), decltype(expr)>::value), "should be true");
            static_assert((IsSubTree<decltype(mv), decltype(expr)>::value), "should be true");
            static_assert((IsSubTree<Matrix, decltype(expr)>::value), "should be true");
            static_assert((IsSubTree<Vector, decltype(expr)>::value), "should be true");
        }

        void trace_test() {
            SizeType n = 3;
            Matrix m;
            m.dense_identity(serial_layout(n, n), 0.5);
            Scalar t = trace(m);
            utopia_test_assert(approxeq(t, 1.5, 1e-16));
        }

        void in_place_test() {
            SizeType n = 3;

            Matrix oracle, m;
            auto mat_layout = layout(world, Traits::decide(), Traits::decide(), n, n);
            oracle.dense_identity(mat_layout);
            m.dense_identity(mat_layout, 0.5);

            // Matrix m = 0.5 * dense_identity(n, n);
            m *= 2.0;

            utopia_test_assert(approxeq(m, oracle, 1e-16));

            oracle *= 1. / 4.0;
            m /= 4.0;

            utopia_test_assert(approxeq(m, oracle, 1e-16));
        }

        void chop_test() {
            SizeType n = 10;
            auto vec_layout = layout(world, Traits::decide(), n);

            Vector x(vec_layout, 1.0);
            Vector y(vec_layout, -1.0);

            {
                auto y_view = view_device(y);
                parallel_for(
                    range_device(y), UTOPIA_LAMBDA(const SizeType &i) {
                        const Scalar val = (i == 0) ? 1e-14 : ((i < n / 2.0) ? -i : i);
                        y_view.set(i, val);
                    });
            }

            Matrix M = outer(x, y);
            chop_smaller_than(M, 1e-13);
            chop_greater_than(M, 1e-13);

            utopia_test_assert(approxeq(0.0, norm2(M), 1e-13));
        }

    public:
        void run() {
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
            UTOPIA_RUN_TEST(chop_test);
        }
    };

    static void algebra() {
#ifdef WITH_BLAS
        DenseAlgebraTest<BlasMatrixd, BlasVectord>().run();
        SerialAlgebraTest<BlasMatrixd, BlasVectord>().run();
        SparseAlgebraTest<BlasMatrixd, BlasVectord>().run();
#endif  // WITH_BLAS

#ifdef UTOPIA_WITH_PETSC
        DenseAlgebraTest<PetscMatrix, PetscVector>().run();
        SparseAlgebraTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC

#ifdef WITH_TRILINOS
        VectorAlgebraTest<TpetraVector>().run();
        SparseAlgebraTest<TpetraMatrix, TpetraVector>().run();
#endif  // WITH_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(algebra);

}  // namespace utopia

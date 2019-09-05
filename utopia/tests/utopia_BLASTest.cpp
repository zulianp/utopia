/*
 * @Author: Eric Botter
 * @Date:   2016-11-15
 */
#include "utopia.hpp"
#include "utopia_BLASTest.hpp"
#include "test_problems/utopia_TestFunctions2D.hpp"
#include "test_problems/utopia_TestFunctionsND.hpp"

typedef double Real;

namespace utopia {

#ifdef WITH_BLAS

    Matrixd hm_matrix(const SizeType rows, const SizeType cols, const std::vector<Real> &values)
    {
        Matrixd mat = zeros(rows, cols);
        mat.implementation().entries() = values;
        return mat;
    }

    Vectord hm_vector(const std::vector<Real> &values)
    {
        Vectord vec = zeros(values.size());
        vec.implementation() = values;
        return vec;
    }

    void blas_gemm_test()
    {
        Matrixd A{ hm_matrix(2, 2, {1, 1, 1, 1}) };
        Matrixd B{ hm_matrix(2, 2, {2, 2, 2, 2}) };
        Matrixd C{ hm_matrix(2, 2, {4, 0, 0, 4}) };
        Matrixd res;

        //specialization 1
        res = 0.1 * A * B;
        Matrixd expected = values(2, 2, 0.4);

        utopia_test_assert(approxeq(expected, res));

        res *= -.6;

        //specialization 2
        res = 0.1 * A * B + 0.2 * C;

        {
            Write<Matrixd> w_r(expected);
            expected.set(0, 0, 1.2);
            expected.set(1, 1, 1.2);
        }

        utopia_test_assert(approxeq(expected, res));
    }

    void blas_test() {
        //variables
        Matrixd m1{ hm_matrix(2, 2, {1, 0, 0, 1}) };
        Matrixd m2{ hm_matrix(2, 2, {1, 2, 3, 4}) };
        Vectord v{ hm_vector({1.0, 10.0}) };

        { //BLAS 2 + 3
            Vectord vresult;
            auto expr = m2 * v - (m1 * m2) * v;
            // std::cout << tree_format(expr.getClass()) << std::endl;
            vresult = expr;

            Vectord vexp({0.0, 0.0});
            utopia_test_assert(approxeq(vexp, vresult));
        }

        Matrixd mexp{ hm_matrix(2, 2, {1, 3, 2, 4}) };

        { //BLAS 3
            Matrixd mresult;
            auto mexpr = transpose(m2) * transpose(m1);
            // std::cout << tree_format(mexpr.getClass()) << std::endl;
            mresult = mexpr;

            utopia_test_assert(approxeq(mexp, mresult));
        }

        {
            Matrixd mresult;
            mresult = transpose(m2);
            utopia_test_assert(approxeq(mexp, mresult));
        }
    }


    void blas_pow_test() {
        //variables

        Vectord v{ hm_vector({2.0, 10.0}) };

        // int a=8;

        utopia::Vectord res=power(v,4.0);

        //auto res2=powAB(2.0);

        //disp(res2);
    }




    void blas_function_test() {
        Vectord point({1.0, -1.0});

        TestFunction2D_1<Matrixd, Vectord> fun;
        Vectord g;
        Matrixd H;

        fun.gradient(point, g);
        fun.hessian(point, H);

        Vectord::Scalar fun_value;
        fun.value(point, fun_value);

        Vectord::Scalar val_exp = 169.0;
        utopia_test_assert(approxeq(val_exp, fun_value));

        Vectord g_exp{ hm_vector({-10.0, 48.0}) };
        utopia_test_assert(approxeq(g_exp, g));

        Matrixd H_exp{ hm_matrix(2, 2, {4.0, 0.0, 0.0, 8.0}) };
        utopia_test_assert(approxeq(H_exp, H));
    }

    void blas_solver_test() {
#ifdef WITH_LAPACK
        auto lapackSolver = std::make_shared< LUDecomposition<Matrixd, Vectord> >();
        Vectord x0({3.0, -2.0});

        Newton<Matrixd, Vectord> newtonSolver(lapackSolver);
        newtonSolver.enable_differentiation_control(false);

        TestFunctionND_1<Matrixd, Vectord> fun2(10);

        x0 = values(10, 2.0);
        newtonSolver.solve(fun2, x0);
#endif //WITH_LAPACK

    }

    void blas_inplace_test() {
        Vectord v1{ hm_vector({4.0, 3.0, 2.0, 1.0}) };
        Vectord v2{ hm_vector({1.0, 2.0, 3.0, 4.0}) };
        Matrixd m1{ hm_matrix(2, 2, {2.0, 2.0, 2.0, 2.0}) };
        // Matrixd m2(2, 2, {2.0, 1.0, 1.0, 1.0});

        //! [in place operations (blas)]
        v1 -= v2;
        m1 *= m1;
        //! [in place operations (blas)]

        Vectord v_exp{ hm_vector({3.0, 1.0, -1.0, -3.0}) };
        utopia_test_assert(approxeq(v_exp, v1));

        Matrixd m_exp{ hm_matrix(2, 2, {8.0, 8.0, 8.0, 8.0}) };
        utopia_test_assert(approxeq(m_exp, m1));
    }

    void blas_accessors_test() {
        Vectord v1{ hm_vector({0.0, 0.0}) };
        Vectord v2{ hm_vector({1.0, 1.0}) };
        Matrixd m1{ hm_matrix(2, 2, {2.0, 2.0, 2.0, 2.0}) };
        Matrixd m2{ hm_matrix(2, 2, {0.0, 0.0, 0.0, 0.0}) };

        {
            Write<Vectord> w_v1(v1);
            v1.set(0, 1);
            v1.set(1, 2);
        }

        {
            Write<Matrixd> w_m2(m2);
            Read<Vectord> r_v1(v1);

            m2.set(0, 0, v1.get(0));
            m2.set(1, 1, 1);
        }

        v1 -= v2;
        m1 *= m2;

        Vectord v_exp{ hm_vector({0.0, 1.0}) };
        utopia_test_assert(approxeq(v_exp, v1));

        Matrixd m_exp{ hm_matrix(2, 2, {2.0, 2.0, 2.0, 2.0}) };
        utopia_test_assert(approxeq(m_exp, m1));



        v1.set(6.);
    }


    void blas_set_values_test() {
        Matrixd m1 = identity(3, 3);

        std::vector<SizeType> rows{0, 1, 2};
        std::vector<SizeType> cols{0, 1, 2};
        std::vector<double> vals{
            100, 0, 0,
            0, 101, 0,
            0, 0, 102};

        {
            Write<Matrixd> w_m1(m1);
            m1.set_matrix(rows, cols, vals);
        }

        Matrixd m_exp{ hm_matrix(3, 3, {100.0, 0, 0, 0, 101.0, 0, 0, 0, 102.0}) };
        utopia_test_assert(approxeq(m_exp, m1));
    }

    void blas_axpy_test()
    {
        Vectord w1{ hm_vector({1., 2., 3.}) };
        Vectord w2{ hm_vector({4., 5., 6.}) };
        Vectord actual = w1 + 0.1 * w2;
        Vectord expected{ hm_vector({1.4, 2.5, 3.6}) };

        utopia_test_assert(approxeq(expected, actual));
    }

    void blas_norm_test() {
        Vectord w1{ hm_vector({10.0, 3.0, 1.0}) };
        Vectord w2{ hm_vector({20.0, 2.0, 3.0}) };
        Vectord w3{ hm_vector({-30.0, -5.0, -4.0}) };
        Vectord wresult;

        auto twiceaxpy = 0.9 * (w1 * 0.1 + w2) + w3; //axpy twice

        wresult = twiceaxpy;
        Vectord wexp{ hm_vector({-11.1, -2.93, -1.21}) };

        utopia_test_assert(approxeq(wexp, wresult));

        Real val = norm2(twiceaxpy);
        // std::cout << tree_format(norm2(twiceaxpy).getClass()) << std::endl;
        // std::cout << val << std::endl;
        val = norm_infty(twiceaxpy); UTOPIA_UNUSED(val);
    }

    void blas_composite_test() {
        Vectord w1{ hm_vector({10.0, 3.0, 1.0}) };
        Vectord w2{ hm_vector({20.0, 2.0, 3.0}) };
        Vectord w3{ hm_vector({-30.0, -5.0, -4.0}) };
        Vectord wresult;
        //advanced (To make it work for all backends)
        auto twiceaxpy = 0.9 * (w1 * 0.1 + w2) + w3; //axpy twice
        auto temp = twiceaxpy + w1 * dot(w1, 4.0 * (2.0 * w2) + w3 * 6.0);
        //    //multiline calculations
        auto expr = temp * 0.01;

        //query the expression structure
        // std::cout << tree_format(expr.getClass()) << std::endl;

        //Evaluate and verify value of the expression
        wresult = expr;
        Vectord wexp{ hm_vector({-24.311, -7.2893, -2.4321}) };
        utopia_test_assert(approxeq(wexp, wresult));
    }

    void blas_row_view_test()
    {
        CRSMatrixd mat = sparse(3, 3, 1);

        {
            Write<CRSMatrixd> write(mat);
            mat.set(0, 0, 2);
            mat.set(1, 1, 2);
            mat.set(2, 2, 2);
        }

        //controlled way
        {
            Read<CRSMatrixd> r_m(mat);

            Size s = size(mat);
            Range r = row_range(mat);
            for(auto i = r.begin(); i != r.end(); ++i) {
                RowView<const CRSMatrixd> row_view(mat, i);
                const SizeType n_values = row_view.n_values();

                utopia_test_assert(n_values == SizeType(1));

                for(SizeType k = 0; k < n_values; ++k) {
                    auto c = row_view.col(k);
                    auto v = row_view.get(k);

                    utopia_test_assert(approxeq(2., v));
                    utopia_test_assert(SizeType(c) == SizeType(i));
                }
            }
        }

        //simple way
        each_read(mat, [](const SizeType i, const SizeType j, const double v) {
            utopia_test_assert(i == j);
            utopia_test_assert(approxeq(2., v));
        });

        Matrixd d_mat = values(3, 3, 2.);

        SizeType n_vals = 0;
        each_read(d_mat, [&n_vals](const SizeType /*i*/, const SizeType /*j*/, const double v) {
            utopia_test_assert(approxeq(2., v));
            ++n_vals;
        });

        utopia_test_assert(n_vals == 3 * 3);
    }


    template<typename SizeType, typename Scalar>
    utopia::CRSMatrixd crs(const SizeType rows, const SizeType cols,
                           std::initializer_list <SizeType> rowptr,
                           std::initializer_list <SizeType> colindex,
                           std::initializer_list <Scalar> values
                           ) {
        utopia::CRSMatrixd ret;
        ret.implementation().initialize(rows, cols, rowptr, colindex, values);
        return ret;
    }

    void blas_sparse_matrix_test() {
        CRSMatrixd mat = sparse(3, 3, 1);

        {
            Write<CRSMatrixd> write(mat);
            mat.set(0, 0, 2);
            mat.set(1, 1, 2);
            mat.set(2, 2, 2);
        }

        Vectord v1{ hm_vector({2, 2, 2}) };
        Vectord v2 = mat * v1;
        utopia_test_assert(approxeq(2 * v1, v2));

        CCSMatrixd ccsmat = mat;

#ifdef WITH_UMFPACK
        //     auto lsolver = std::make_shared<UmfpackLU>();

        //     lsolver->solve(mat, v1, v2);

        //     Newton<CRSMatrixd, Vectord> nlsolver(lsolver);
        // //    nlsolver.enable_differentiation_control(true);

        //     SimpleQuadraticFunction<CRSMatrixd, Vectord> fun;
        //     Vectord x = values(100, 2.);
        //     Vectord expected = zeros(x.size());

        //     nlsolver.solve(fun, x);
        //     utopia_test_assert(approxeq(expected, x));


        //     Rosenbrock01<CRSMatrixd, Vectord> rosenbrock;
        //     Vectord x2({-2, -2});
        //     nlsolver.solve(rosenbrock, x2);

        //     expected = values(2, 1);
        //     utopia_test_assert(approxeq(expected, x2));
#endif //WITH_UMFPACK
    }

    void blas_pgs_test()
    {
        int n = 3;
        CRSMatrixd A = sparse(n, n, 1);

        {
            Write<CRSMatrixd> w_A(A);
            A.set(0, 0, 1);
            A.set(1, 1, 1);
            A.set(2, 2, 1);
        }

        Vectord rhs = values(n, 2.);
        Vectord x   = zeros(n);

        ProjectedGaussSeidel<CRSMatrixd, Vectord> pgs;
        pgs.solve(A, rhs, x);
    }

#endif //WITH_BLAS

    void runBLASTest() {
#ifdef WITH_BLAS
        UTOPIA_UNIT_TEST_BEGIN("BLASTest");
        UTOPIA_RUN_TEST(blas_pgs_test);
        UTOPIA_RUN_TEST(blas_gemm_test);
        UTOPIA_RUN_TEST(blas_row_view_test);
        UTOPIA_RUN_TEST(blas_test);
        UTOPIA_RUN_TEST(blas_axpy_test);
        UTOPIA_RUN_TEST(blas_function_test);
        UTOPIA_RUN_TEST(blas_solver_test);
        UTOPIA_RUN_TEST(blas_inplace_test);
        UTOPIA_RUN_TEST(blas_accessors_test);
        UTOPIA_RUN_TEST(blas_set_values_test);
        UTOPIA_RUN_TEST(blas_norm_test);
        UTOPIA_RUN_TEST(blas_composite_test);
        UTOPIA_RUN_TEST(blas_sparse_matrix_test);
        UTOPIA_RUN_TEST(blas_pow_test);
        UTOPIA_UNIT_TEST_END("BLASTest");
#endif // WITH_BLAS
    }
}

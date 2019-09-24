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

    BlasMatrixd hm_matrix(const SizeType rows, const SizeType cols, const std::vector<Real> &values)
    {
        BlasMatrixd mat;
        mat.resize(rows, cols);
        mat.entries() = values;
        return mat;
    }

    BlasVectord hm_vector(const std::vector<Real> &values)
    {
        BlasVectord vec;
        vec.resize(values.size());
        vec.entries() = values;
        return vec;
    }

    void blas_gemm_test()
    {
        BlasMatrixd A{ hm_matrix(2, 2, {1, 1, 1, 1}) };
        BlasMatrixd B{ hm_matrix(2, 2, {2, 2, 2, 2}) };
        BlasMatrixd C{ hm_matrix(2, 2, {4, 0, 0, 4}) };
        BlasMatrixd res;

        //specialization 1
        res = 0.1 * A * B;
        BlasMatrixd expected = values(2, 2, 0.4);

        utopia_test_assert(approxeq(expected, res));

        res *= -.6;

        //specialization 2
        res = 0.1 * A * B + 0.2 * C;

        {
            Write<BlasMatrixd> w_r(expected);
            expected.set(0, 0, 1.2);
            expected.set(1, 1, 1.2);
        }

        utopia_test_assert(approxeq(expected, res));
    }

    void blas_test() {
        //variables
        BlasMatrixd m1{ hm_matrix(2, 2, {1, 0, 0, 1}) };
        BlasMatrixd m2{ hm_matrix(2, 2, {1, 2, 3, 4}) };
        BlasVectord v{ hm_vector({1.0, 10.0}) };

        { //BLAS 2 + 3
            BlasVectord vresult;
            auto expr = m2 * v - (m1 * m2) * v;
            // std::cout << tree_format(expr.get_class()) << std::endl;
            vresult = expr;

            BlasVectord vexp({0.0, 0.0});
            utopia_test_assert(approxeq(vexp, vresult));
        }

        BlasMatrixd mexp{ hm_matrix(2, 2, {1, 3, 2, 4}) };

        { //BLAS 3
            BlasMatrixd mresult;
            auto mexpr = transpose(m2) * transpose(m1);
            // std::cout << tree_format(mexpr.get_class()) << std::endl;
            mresult = mexpr;

            utopia_test_assert(approxeq(mexp, mresult));
        }

        {
            BlasMatrixd mresult;
            mresult = transpose(m2);
            utopia_test_assert(approxeq(mexp, mresult));
        }
    }


    void blas_pow_test() {
        //variables

        BlasVectord v{ hm_vector({2.0, 10.0}) };

        // int a=8;

        utopia::BlasVectord res=power(v,4.0);

        //auto res2=powAB(2.0);

        //disp(res2);
    }

    void blas_function_test() {
        BlasVectord point({1.0, -1.0});

        TestFunction2D_1<BlasMatrixd, BlasVectord> fun;
        BlasVectord g;
        BlasMatrixd H;

        fun.gradient(point, g);
        fun.hessian(point, H);

        BlasVectord::Scalar fun_value;
        fun.value(point, fun_value);

        BlasVectord::Scalar val_exp = 169.0;
        utopia_test_assert(approxeq(val_exp, fun_value));

        BlasVectord g_exp{ hm_vector({-10.0, 48.0}) };
        utopia_test_assert(approxeq(g_exp, g));

        BlasMatrixd H_exp{ hm_matrix(2, 2, {4.0, 0.0, 0.0, 8.0}) };
        utopia_test_assert(approxeq(H_exp, H));
    }

    void blas_solver_test() {
#ifdef WITH_LAPACK
        auto lapackSolver = std::make_shared< LUDecomposition<BlasMatrixd, BlasVectord> >();
        BlasVectord x0({3.0, -2.0});

        Newton<BlasMatrixd, BlasVectord> newtonSolver(lapackSolver);
        newtonSolver.enable_differentiation_control(false);

        TestFunctionND_1<BlasMatrixd, BlasVectord> fun2(10);

        x0 = values(10, 2.0);
        newtonSolver.solve(fun2, x0);
#endif //WITH_LAPACK

    }

    void blas_inplace_test() {
        BlasVectord v1{ hm_vector({4.0, 3.0, 2.0, 1.0}) };
        BlasVectord v2{ hm_vector({1.0, 2.0, 3.0, 4.0}) };
        BlasMatrixd m1{ hm_matrix(2, 2, {2.0, 2.0, 2.0, 2.0}) };
        // BlasMatrixd m2(2, 2, {2.0, 1.0, 1.0, 1.0});

        //! [in place operations (blas)]
        v1 -= v2;
        m1 *= m1;
        //! [in place operations (blas)]

        BlasVectord v_exp{ hm_vector({3.0, 1.0, -1.0, -3.0}) };
        utopia_test_assert(approxeq(v_exp, v1));

        BlasMatrixd m_exp{ hm_matrix(2, 2, {8.0, 8.0, 8.0, 8.0}) };
        utopia_test_assert(approxeq(m_exp, m1));
    }

    void blas_accessors_test() {
        BlasVectord v1{ hm_vector({0.0, 0.0}) };
        BlasVectord v2{ hm_vector({1.0, 1.0}) };
        BlasMatrixd m1{ hm_matrix(2, 2, {2.0, 2.0, 2.0, 2.0}) };
        BlasMatrixd m2{ hm_matrix(2, 2, {0.0, 0.0, 0.0, 0.0}) };

        {
            Write<BlasVectord> w_v1(v1);
            v1.set(0, 1);
            v1.set(1, 2);
        }

        {
            Write<BlasMatrixd> w_m2(m2);
            Read<BlasVectord> r_v1(v1);

            m2.set(0, 0, v1.get(0));
            m2.set(1, 1, 1);
        }

        v1 -= v2;
        m1 *= m2;

        BlasVectord v_exp{ hm_vector({0.0, 1.0}) };
        utopia_test_assert(approxeq(v_exp, v1));

        BlasMatrixd m_exp{ hm_matrix(2, 2, {2.0, 2.0, 2.0, 2.0}) };
        utopia_test_assert(approxeq(m_exp, m1));



        v1.set(6.);
    }


    void blas_set_values_test() {
        using SizeType = Traits<BlasMatrixd>::SizeType;

        BlasMatrixd m1 = identity(3, 3);

        std::vector<SizeType> rows{0, 1, 2};
        std::vector<SizeType> cols{0, 1, 2};
        std::vector<double> vals{
            100, 0, 0,
            0, 101, 0,
            0, 0, 102};

        {
            Write<BlasMatrixd> w_m1(m1);
            m1.set_matrix(rows, cols, vals);
        }

        BlasMatrixd m_exp{ hm_matrix(3, 3, {100.0, 0, 0, 0, 101.0, 0, 0, 0, 102.0}) };
        utopia_test_assert(approxeq(m_exp, m1));
    }

    void blas_axpy_test()
    {
        BlasVectord w1{ hm_vector({1., 2., 3.}) };
        BlasVectord w2{ hm_vector({4., 5., 6.}) };
        w1 += 0.1 * w2;
        BlasVectord expected{ hm_vector({1.4, 2.5, 3.6}) };

        utopia_test_assert(approxeq(expected, w1));
    }

    void blas_norm_test() {
        BlasVectord w1{ hm_vector({10.0, 3.0, 1.0}) };
        BlasVectord w2{ hm_vector({20.0, 2.0, 3.0}) };
        BlasVectord w3{ hm_vector({-30.0, -5.0, -4.0}) };
        BlasVectord wresult;

        auto twiceaxpy = 0.9 * (w1 * 0.1 + w2) + w3; //axpy twice

        wresult = twiceaxpy;
        BlasVectord wexp{ hm_vector({-11.1, -2.93, -1.21}) };

        utopia_test_assert(approxeq(wexp, wresult));

        Real val = norm2(twiceaxpy);
        // std::cout << tree_format(norm2(twiceaxpy).get_class()) << std::endl;
        // std::cout << val << std::endl;
        val = norm_infty(twiceaxpy); UTOPIA_UNUSED(val);
    }

    void blas_composite_test() {
        BlasVectord w1{ hm_vector({10.0, 3.0, 1.0}) };
        BlasVectord w2{ hm_vector({20.0, 2.0, 3.0}) };
        BlasVectord w3{ hm_vector({-30.0, -5.0, -4.0}) };
        BlasVectord wresult;
        //advanced (To make it work for all backends)
        auto twiceaxpy = 0.9 * (w1 * 0.1 + w2) + w3; //axpy twice
        auto temp = twiceaxpy + w1 * dot(w1, 4.0 * (2.0 * w2) + w3 * 6.0);
        //    //multiline calculations
        auto expr = temp * 0.01;

        //query the expression structure
        // std::cout << tree_format(expr.get_class()) << std::endl;

        //Evaluate and verify value of the expression
        wresult = expr;
        BlasVectord wexp{ hm_vector({-24.311, -7.2893, -2.4321}) };
        utopia_test_assert(approxeq(wexp, wresult));
    }

    void blas_row_view_test()
    {
        using MatrixT  = utopia::BlasMatrixd;
        using SizeType = Traits<MatrixT>::SizeType;

        SizeType n = 3;
        MatrixT mat = zeros(n, n);

        {
            Write<MatrixT> write(mat);
            
            for(SizeType i = 0; i < n; ++i) {
                mat.set(i, i, 2);
            }
        }


        //controlled way
        {
            Read<MatrixT> r_m(mat);

            Size s = size(mat);
            Range r = row_range(mat);
            for(auto i = r.begin(); i != r.end(); ++i) {
                RowView<const MatrixT> row_view(mat, i);
                const SizeType n_values = row_view.n_values();

                utopia_test_assert(n_values == n);

                for(SizeType k = 0; k < n_values; ++k) {
                    const SizeType c = row_view.col(k);
                    auto v = row_view.get(k);

                    if(c == SizeType(i)) {
                        utopia_test_assert(approxeq(2., v));
                    }
                }
            }
        }

        //simple way
        each_read(mat, [](const SizeType i, const SizeType j, const double v) {
            if(i==j) {
                utopia_test_assert(approxeq(2., v));
            } else {
                utopia_test_assert(approxeq(0., v));
            }
        });

        BlasMatrixd d_mat = values(3, 3, 2.);

        SizeType n_vals = 0;
        each_read(d_mat, [&n_vals](const SizeType /*i*/, const SizeType /*j*/, const double v) {
            utopia_test_assert(approxeq(2., v));
            ++n_vals;
        });

        utopia_test_assert(n_vals == 3 * 3);
    }

    void blas_pgs_test()
    {
        int n = 3;
        BlasMatrixd A = zeros(n, n);

        {
            Write<BlasMatrixd> w_A(A);
            A.set(0, 0, 1);
            A.set(1, 1, 1);
            A.set(2, 2, 1);
        }

        BlasVectord rhs = values(n, 2.);
        BlasVectord x   = zeros(n);

        ProjectedGaussSeidel<BlasMatrixd, BlasVectord> pgs;
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
        UTOPIA_RUN_TEST(blas_pow_test);
        UTOPIA_UNIT_TEST_END("BLASTest");
#endif // WITH_BLAS
    }
}

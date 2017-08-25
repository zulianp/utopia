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

    void BLAS_test() {
        //variables
        Matrixd m1(2, 2, {1, 0, 0, 1});
        Matrixd m2(2, 2, {1, 2, 3, 4});
        Vectord v({1.0, 10.0});

        { //BLAS 2 + 3
            Vectord vresult;
            auto expr = m2 * v - (m1 * m2) * v;
            // std::cout << tree_format(expr.getClass()) << std::endl;
            vresult = expr;

            Vectord vexp({0.0, 0.0});
            assert(approxeq(vexp, vresult));
        }

        Matrixd mexp(2, 2, {1, 3, 2, 4});

        { //BLAS 3
            Matrixd mresult;
            auto mexpr = transpose(m2) * transpose(m1);
            // std::cout << tree_format(mexpr.getClass()) << std::endl;
            mresult = mexpr;

            assert(approxeq(mexp, mresult));
        }

        {
            Matrixd mresult;
            mresult = transpose(m2);
            assert(approxeq(mexp, mresult));
        }
    }


    void function_test() {
        Vectord point({1.0, -1.0});

        TestFunction2D_1<Matrixd, Vectord> fun;
        Vectord g;
        Matrixd H;

        fun.gradient(point, g);
        fun.hessian(point, H);

        Vectord::Scalar fun_value;
        fun.value(point, fun_value);

        Vectord::Scalar val_exp = 169.0;
        assert(approxeq(val_exp, fun_value));

        Vectord g_exp({-10.0, 48.0});
        assert(approxeq(g_exp, g));

        Matrixd H_exp(2, 2, {4.0, 0.0, 0.0, 8.0});
        assert(approxeq(H_exp, H));
    }

    void solver_test() {
    #ifdef WITH_LAPACK
    //    std::cout << "begin: solver_test" << std::endl;
    //    Matrixd A(5, 5,
    //              {6.80, -2.11, 5.66, 5.97, 8.23,
    //               -6.05, -3.30, 5.36, -4.44, 1.08,
    //               -0.45, 2.58, -2.70, 0.27, 9.04,
    //               8.32, 2.71, 4.35, -7.17, 2.14,
    //               -9.67, -5.14, -7.26, 6.08, -6.87});
    //    Vectord b({4.02, 6.19, -8.22, -7.57, -3.03});
    //    Vectord x({0, 0, 0, 0, 0});
    //
        auto lapackSolver = std::make_shared< LUDecomposition<Matrixd, Vectord> >();
    //
    //    lapackSolver->solve(A.implementation(), b.implementation(), x.implementation());
    //    disp(x);
    //
    //    const SimpleQuadraticFunction<Matrixd, Vectord> fun;
        Vectord x0({3.0, -2.0});
    //
        Newton<Matrixd, Vectord> newtonSolver(lapackSolver);
        newtonSolver.enable_differentiation_control(false);

    //
    //    disp(x0);
    //
        // TrustRegion<Matrixd, Vectord> trustRegionSolver(lapackSolver);
        // trustRegionSolver.enable_differentiation_control(false);
    //    Vectord x1({-1.0, -2.0});
    //
        TestFunctionND_1<Matrixd, Vectord> fun2(10);
    //
        x0 = values(10, 2.0);
    //    APTS_2Domains<Matrixd, Vectord> APTSSolver(lapackSolver);
    //    APTSSolver.solve(fun2, x0);

        newtonSolver.solve(fun2, x0);

    //
    //    disp(x0);
    //    std::cout << fun2.value(x0) << std::endl;
    //
    //    const Rosenbrock<Matrixd, Vectord> rosenbrock;
    //    trustRegionSolver.solve(rosenbrock, x1);
    //    disp(x1);

    //    x0 = values(10, 2.0);
    //    trustRegionSolver.solve(fun2, x0);
    //    disp(x0);


    //    const unsigned FE_nodes = 7;
    //
    //    NeoHookean1D<Matrixd, Vectord> fun3(FE_nodes, 1, 1);
    //
    //    Vectord x3;
    //    x3 = values(FE_nodes, 0.0);
    //    trustRegionSolver.solve(fun3, x3);
    //
    //    disp(x3);
    //    disp(fun3.value(x3));
    //
    //    APTS2<Matrixd, Vectord> APTSSolver2(lapackSolver);
    //    x3 = values(FE_nodes, 0.0);
    //    x3.set(0, 0.0);
    //    APTSSolver2.solve(fun3, x3);
    //
    //    disp(x3);
    //    disp(fun3.value(x3));
        // std::cout << "end: solver_test" << std::endl;
    #endif //WITH_LAPACK
    }

    void inplace_test() {
        // std::cout << "begin: inplace_test" << std::endl;
        //! [in place operations (blas)]

        Vectord v1({4.0, 3.0, 2.0, 1.0});
        Vectord v2({1.0, 2.0, 3.0, 4.0});
        Matrixd m1(2, 2, {2.0, 2.0, 2.0, 2.0});
        // Matrixd m2(2, 2, {2.0, 1.0, 1.0, 1.0});

        v1 -= v2;
        m1 *= m1;

        Vectord v_exp({3.0, 1.0, -1.0, -3.0});
        assert(approxeq(v_exp, v1));

        Matrixd m_exp(2, 2, {8.0, 8.0, 8.0, 8.0});
        assert(approxeq(m_exp, m1));

        //! [in place operations (blas)]
        // std::cout << "end: inplace_test" << std::endl;
    }

    void accessors_test() {
        // std::cout << "begin: accessors_test" << std::endl;

        Vectord v1({0.0, 0.0});
        Vectord v2({1.0, 1.0});
        Matrixd m1(2, 2, {2.0, 2.0, 2.0, 2.0});
        Matrixd m2(2, 2, {0.0, 0.0, 0.0, 0.0});

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

        Vectord v_exp({0.0, 1.0});
        assert(approxeq(v_exp, v1));

        Matrixd m_exp(2, 2, {2.0, 2.0, 2.0, 2.0});
        assert(approxeq(m_exp, m1));

        // std::cout << "end: accessors_test" << std::endl;
    }


    void setvalues_test() {
        // std::cout << "begin: setvalues_test" << std::endl;
        Matrixd m1 = identity(3, 3);

        std::vector<int> rows{0, 1, 2};
        std::vector<int> cols{0, 1, 2};
        std::vector<double> vals{100, 101, 102};

        {
            Write<Matrixd> w_m1(m1);
            m1.set(rows, cols, vals);
        }

        Matrixd m_exp(3, 3, {100.0, 0, 0, 0, 101.0, 0, 0, 0, 102.0});
        assert(approxeq(m_exp, m1));

        // std::cout << "end: setvalues_test" << std::endl;
    }

    void norm_test() {
        // std::cout << "begin: norm_test" << std::endl;

        Vectord w1({10.0, 3.0, 1.0});
        Vectord w2({20.0, 2.0, 3.0});
        Vectord w3({-30.0, -5.0, -4.0});
        Vectord wresult;

        auto twiceaxpy = 0.9 * (w1 * 0.1 + w2) + w3; //axpy twice

        wresult = twiceaxpy;
        Vectord wexp({-11.1, -2.93, -1.21});
        assert(approxeq(wexp, wresult));

        Real val = norm2(twiceaxpy);
        // std::cout << tree_format(norm2(twiceaxpy).getClass()) << std::endl;
        // std::cout << val << std::endl;
        val = norm_infty(twiceaxpy);

        // std::cout << "end: norm_test" << std::endl;
    }

    void composite_test() {
        // std::cout << "begin: composite_test" << std::endl;

        Vectord w1({10.0, 3.0, 1.0});
        Vectord w2({20.0, 2.0, 3.0});
        Vectord w3({-30.0, -5.0, -4.0});
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
        Vectord wexp({-24.311, -7.2893, -2.4321});
        assert(approxeq(wexp, wresult));

        // std::cout << "end: composite_test" << std::endl;
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

    void sparse_matrix_test() {
        CRSMatrixd mat = sparse(3, 3, 1);

        {
            Write<CRSMatrixd> write(mat);
            mat.set({0, 1, 2}, {0, 1, 2}, {2, 2, 2});
        }

        Vectord v1({2, 2, 2});
        Vectord v2 = mat * v1;
        assert(approxeq(2 * v1, v2));

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
        //     assert(approxeq(expected, x));


        //     Rosenbrock<CRSMatrixd, Vectord> rosenbrock;
        //     Vectord x2({-2, -2});
        //     nlsolver.solve(rosenbrock, x2);

        //     expected = values(2, 1);
        //     assert(approxeq(expected, x2));
        #endif //WITH_UMFPACK

    }

    #endif //WITH_BLAS

    void runBLASTest() {
        #ifdef WITH_BLAS
        std::cout << "Begin: BLASTest" << std::endl;

        BLAS_test();
        function_test();
        solver_test();
        inplace_test();
        accessors_test();
        setvalues_test();
        norm_test();
        composite_test();
        sparse_matrix_test();

        std::cout << "End:   BLASTest" << std::endl;
        #endif // WITH_BLAS
    }
}

/*
* @Author: alenakopanicakova
* @Date:   2016-07-15
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-02-07
*/
#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"

namespace utopia
{
     /**
      * @brief      Class to test our nonlinear solvers.
      *
      * @todo       add more tests to TR, LS, ...
      * @todo       different strategies, different LS solve, test all params, nonlinear problems ...
      *
      * @tparam     Matrix
      * @tparam     Vector
      * @tparam     Scalar
      */
    template<class Matrix, class Vector, class Scalar>
    class SolverTest
    {
     public:

        void run()
        {
            newton_cg_utopia_test();
            solver_from_params();
            TR_test();
            LS_test();
            nl_solve_test();
        }


        class EmptyLSFun : public LeastSquaresFunction<Matrix, Vector> {
        public:
            bool residual(const Vector &/*point*/, Vector &/*result*/) const {}
            bool jacobian(const Vector &/*x*/, Matrix &/*hessian*/) const {}
            bool value(const Vector &, Scalar &val) const {}
            bool update(const Vector &) { }
        };

        void ls_normal_eq()
        {
            LeastSquaresNewton<Matrix, Vector> newton(std::make_shared<ConjugateGradient<Matrix, Vector>>());
            auto ls_strat  = std::make_shared<utopia::Backtracking<Matrix, Vector> >();
            newton.set_line_search_strategy(ls_strat);
            
            EmptyLSFun fun;
            Vector x0;
            newton.solve(fun, x0);
        }

        void nl_solve_test()
        {
            //! [NL solve example]
            using namespace utopia;

            //set-up problem
            int n = 10;
            Vector actual = values(n, 2.);
            TestFunctionND_1<Matrix, Vector> fun(n);

            //solve problem
            solve(fun, actual);

            //test outcome...
            Vector expected = values(n, 0.468919);
            assert(approxeq(expected, actual));
            //! [NL solve example]
        }

        void newton_cg_utopia_test()
        {
            //! [Newton CG example]
            using namespace utopia;
            using namespace std;

            //CG with diagonal preconditioner
            auto linear_solver  = make_shared< ConjugateGradient<Matrix, Vector> >();
            auto preconditioner = make_shared< InvDiagPreconditioner<Matrix, Vector> >();
            linear_solver->set_preconditioner(preconditioner);

            //Newton solver with cg linear solver
            Newton<Matrix, Vector> newton_solver(linear_solver);

            const int n = 10;
            Vector actual   = values(n, 2.);
            Vector expected = values(n, 0.468919);

            TestFunctionND_1<Matrix, Vector> fun(n);

            newton_solver.solve(fun, actual);

            //Check if the result is what we expected
            assert(approxeq(expected, actual));
            //! [Newton CG example]
        }

        void solver_from_params()
        {
            using namespace utopia;
           // std::cout << "         Begin: solver_from_params" << std::endl;
            Vector x = values(10, 2.0);
            Vector expected = values(x.size().get(0), 0.468919);
            TestFunctionND_1<Matrix, Vector> fun2(x.size().get(0));

            Parameters params;
            params.tol(1e-7);
            params.solver_type(TRUST_REGION_TAG);
            params.lin_solver_type(BICGSTAB_TAG);
            params.linear_solver_verbose(false);
            params.verbose(false);

            // solve(fun2, x, params);

            // monitor test
            // Matrix H = identity(2,2);
            // int i = 0;
            // monitor(i, H);
            // Matrixd blas_mat  = identity(3,3);
            // monitor(i, blas_mat);

           // std::cout << "         End: solver_from_params" << std::endl;
        }

        void TR_test()
        {
            // rosenbrock test
            if(mpi_world_size() == 1)
            {
           // std::cout << "         Begin: TR_test" << std::endl;
            Vector x = values(10, 2);
            TestFunctionND_1<Matrix, Vector> fun2(x.size().get(0));
            Vector expected = values(x.size().get(0), 0.468919);

            Vector x_w1  = values(4, 10);
            Vector expected_woods = values(4, 1);
            Woods<Matrix, Vector> fun_woods;
            {
                Write<Vector> w1(x_w1);
                x_w1.set(0, -3);
                x_w1.set(1, -1);
                x_w1.set(2, -3);
                x_w1.set(3, -1);
            }


            Parameters params;
            params.atol(1e-10);
            params.rtol(1e-10);
            params.stol(1e-10);
            params.solver_type(TRUST_REGION_TAG);
            params.lin_solver_type(BICGSTAB_TAG);
            params.trust_region_alg(DOGLEG_TAG);
            params.verbose(false);

            // trust_region_solve(fun2, x, params);
            // trust_region_solve(fun_woods, x_w1, params);
            // assert(approxeq(expected, x));


            x = values(10, 2);
            params.trust_region_alg(STEIHAUG_TOINT_TAG);
            trust_region_solve(fun2, x, params);
            assert(approxeq(expected, x));


            // x = values(10, 2);
            // params.trust_region_alg(STEIHAUG_TOINT_TAG);

            // auto strategy_tr = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
            // auto precond = std::make_shared< InvDiagPreconditioner<Matrix, Vector> >();
            // strategy_tr->set_preconditioner(precond);
            // strategy_tr->verbose(true);
            // TrustRegion<Matrix, Vector> tr_solver(strategy_tr);

            // tr_solver.verbose(true);
            // tr_solver.solve(fun2, x);
            // assert(approxeq(expected, x));



            x = values(10, 2);
            params.trust_region_alg(CAUCHYPOINT_TAG);
            trust_region_solve(fun2, x, params);
            assert(approxeq(expected, x));

            x = values(10, 2);
            trust_region_solve(fun2, x, params);
            assert(approxeq(expected, x));


                Vector expected_rosenbrock = values(2, 1);

                Rosenbrock<Matrix, Vector> rosenbrock;
                // params.trust_region_alg(DOGLEG_TAG);
                Vector x0 = values(2, 2.0);
                // trust_region_solve(rosenbrock, x0, params);
                // assert(approxeq(expected_rosenbrock, x0));

                x0 = values(2, 2.0);
                params.trust_region_alg(STEIHAUG_TOINT_TAG);
                trust_region_solve(rosenbrock, x0, params);
                assert(approxeq(expected_rosenbrock, x0));



            }

           // std::cout << "         End: TR_test" << std::endl;

        }

        void LS_test()
        {

            if(mpi_world_size() == 1) {
           // std::cout << "         Begin: LS_test" << std::endl;
            Vector x1 = values(10, 2);
            Vector x2 = values(10, 2);
            TestFunctionND_1<Matrix, Vector> fun2(x1.size().get(0));

            Vector expected = values(x1.size().get(0), 0.468919);

            Parameters params ;
            params.atol(1e-11);
            params.rtol(1e-11);
            params.stol(1e-11);
            params.verbose(false);
            params.linear_solver_verbose(false);
            params.line_search_inner_verbose(false);

            auto lsolver = std::make_shared< ConjugateGradient<Matrix, Vector, -1> >();
            Newton<Matrix, Vector> nlsolver1(lsolver);
            Newton<Matrix, Vector> nlsolver2(lsolver);


            auto strategy_sbc = std::make_shared<utopia::SimpleBacktracking<Matrix, Vector> >();
            auto strategy_bc  = std::make_shared<utopia::Backtracking<Matrix, Vector> >();

            strategy_sbc->set_parameters(params);
            strategy_bc->set_parameters(params);

            nlsolver1.set_line_search_strategy(strategy_sbc);
            nlsolver2.set_line_search_strategy(strategy_bc);

            nlsolver1.set_parameters(params);
            nlsolver2.set_parameters(params);

            nlsolver1.solve(fun2, x1);
            nlsolver2.solve(fun2, x2);

            // Woods function test
            Vector x_w1  = values(4, 10);
            Vector x_w2  = values(4, 10);
            Vector expected_woods = values(4, 1);
            {
                Write<Vector> w1(x_w1);
                Write<Vector> w2(x_w2);
                x_w1.set(0, -3);    x_w2.set(0, -3);
                x_w1.set(1, -1);    x_w2.set(1, -1);
                x_w1.set(2, -3);    x_w2.set(2, -3);
                x_w1.set(3, -1);    x_w2.set(3, -1);
            }

            Woods<Matrix, Vector> fun_woods;
            nlsolver1.solve(fun_woods, x_w1);
            nlsolver2.solve(fun_woods, x_w2);

            assert(approxeq(expected_woods, x_w1));
            assert(approxeq(expected_woods, x_w2));

            // rastrigin function test - convergence to local minimum
            Rastrigin<Matrix, Vector> fun_rastrigin;
            Vector x_r1 = values(2, 1), x_r2 = values(2, 1), expected_rastrigin = values(2, 1);
            {
                Write<Vector> w1(x_r1);
                Write<Vector> w2(x_r2);
                Write<Vector> w3(expected_rastrigin);
                x_r1.set(0, -5.12); x_r2.set(0, -5.12);  expected_rastrigin.set(0, -4.97469);
                x_r1.set(1, 5.12); x_r2.set(1, 5.12);  expected_rastrigin.set(1, 4.97469);
            }

            nlsolver1.solve(fun_rastrigin, x_r2);
            nlsolver2.solve(fun_rastrigin, x_r2);

            // rosenbrock test

                Vector expected_rosenbrock = values(2, 1);
                Rosenbrock<Matrix, Vector> rosenbrock_fun;

                Vector x01 = values(2, 2.0), x02 = values(2, 2.0);
                nlsolver1.solve(rosenbrock_fun, x01);
                nlsolver2.solve(rosenbrock_fun, x02);
                assert(approxeq(expected_rosenbrock, x01));
                assert(approxeq(expected_rosenbrock, x02));

                x01 = values(2, 2.0);
                params.verbose(false);
                params.linear_solver_verbose(false);
                line_search_solve(rosenbrock_fun, x01, params);
            }

           // std::cout << "         End: LS_test" << std::endl;
        }


        SolverTest()
        : _n(10) { }

    private:
        int _n;

    };

#ifdef WITH_PETSC
    class PETScSolverTest {
    public:

        void run()
        {
            petsc_bicgstab_test();
            petsc_gmres_test();
            petsc_newton_test();
            petsc_newton_rosenbrock_test();
            petsc_sparse_semismooth_newton_test();
            petsc_sparse_nonlinear_semismooth_newton_test();
            petsc_direct_solver_newton_test();
            petsc_newton_test_outInfo();
            petsc_sparse_newton_test();
            MG_test();
            CG_MG_test();
            PETSC_CG_MG_test();
            petsc_newton_PETScCG_utopia_test();
            petsc_tr_rr_test();
            petsc_newton_inexact_newton_with_KSP_test(); 

            // neohookean_tr_test();
            // QP_example_MG();
            // TR_in_inf_norm_with_MG();
        }





        // void neohookean_tr_test()
        // {
        //     // std::cout<<"begin: neohook tr test \n";

        //     const unsigned FE_nodes     = 10;
        //     const unsigned int mu       = 1;
        //     const unsigned int lambda   = 1;
        //     GlobalAndLocalNeoHookean1D<DSMatrixd, DVectord, Matrixd, Vectord> fun_neohook(FE_nodes, mu, lambda);

        //     Parameters params;
        //     params.atol(1e-8);
        //     params.rtol(1e-8);
        //     params.stol(1e-8);
        //     params.solver_type(TRUST_REGION_TAG);
        //     params.trust_region_alg(STEIHAUG_TOINT_TAG);
        //     params.verbose(false);

        //     DVectord x = values(FE_nodes * mpi_world_size(), 0);

        //     auto subproblem = std::make_shared<utopia::KSP_TR<DSMatrixd, DVectord> >();
        //     TrustRegion<DSMatrixd, DVectord> tr_solver(subproblem);
        //     tr_solver.set_parameters(params);
        //     tr_solver.solve(fun_neohook, x);


        //     // std::cout<<"end: neohooke global tr test \n";
        //     // std::cout<<"begin:  local neohook tr test \n";

        //     auto local_subproblem = std::make_shared<utopia::SteihaugToint<Matrixd, Vectord> >();


        //     // auto lapackSolver = std::make_shared< LUDecomposition<Matrixd, Vectord> >();
        //     Local_TR<DSMatrixd, DVectord, Matrixd, Vectord> local_TR(local_subproblem);
        //     local_TR.set_parameters(params);

        //     // APTS_test
        //     APTS_with_overlap<DSMatrixd, DVectord, Matrixd, Vectord> apts_solver(local_TR, subproblem);
        //     apts_solver.set_parameters(params);

        //     x = values(FE_nodes * mpi_world_size(), 0);
        //     apts_solver.solve(fun_neohook, x);

        //     // std::cout<<"end: local neohooke global tr test \n";
        // }



        void petsc_tr_rr_test()
        {
            // rosenbrock test
            if(mpi_world_size() == 1)
            {
                // std::cout << "         Begin: TR_test_KSP" << std::endl;
                DVectord x = values(10, 2);
                TestFunctionND_1<DMatrixd, DVectord> fun2(x.size().get(0));
                DVectord expected = values(x.size().get(0), 0.468919);

                Parameters params;
                params.atol(1e-10);
                params.rtol(1e-10);
                params.stol(1e-10);
                params.verbose(false);

                auto subproblem = std::make_shared<utopia::KSP_TR<DMatrixd, DVectord> >();
                TrustRegion<DMatrixd, DVectord> tr_solver(subproblem);
                tr_solver.set_parameters(params);
                tr_solver.solve(fun2, x);

                params.trust_region_alg(STEIHAUG_TOINT_TAG);
                x = values(10, 2);
                trust_region_solve(fun2, x, params);

                params.trust_region_alg(NASH_TAG);
                x = values(10, 2);
                trust_region_solve(fun2, x, params);

                params.trust_region_alg(LANCZOS_TAG);
                x = values(10, 2);
                trust_region_solve(fun2, x, params);

                params.trust_region_alg(CGNE_TAG);
                x = values(10, 2);
                trust_region_solve(fun2, x, params);

            }
            // std::cout << "         End: TR_test_KSP" << std::endl;
        }



        void petsc_bicgstab_test()
        {
            using namespace utopia;

           // std::cout << "         Begin: BICGSTAB_test" << std::endl;
            DMatrixd mat = identity(_n, _n);
            DVectord rhs = zeros(_n);
            DVectord sol = zeros(_n);

            BiCGStab<DMatrixd, DVectord> bicgs;
            bicgs.solve(mat, rhs, sol);

            DVectord expected = zeros(_n);
            assert(approxeq(expected, sol));
           // std::cout << "         End: BICGSTAB_test" << std::endl;
        }

        void petsc_gmres_test()
        {
            using namespace utopia;

            // std::cout << "         Begin: BICGSTAB_test" << std::endl;
            DMatrixd mat = identity(_n, _n);
            DVectord rhs = zeros(_n);
            DVectord sol = zeros(_n);

            GMRES<DMatrixd, DVectord> gmres;
            gmres.solve(mat, rhs, sol);

            DVectord expected = zeros(_n);
            assert(approxeq(expected, sol));
           // std::cout << "         End: BICGSTAB_test" << std::endl;
        }

        void petsc_newton_test_outInfo()
        {
            using namespace utopia;

            // std::cout << "         Begin: petsc_newton_test_outInfo" << std::endl;
            auto lsolver = std::make_shared< BiCGStab<DMatrixd, DVectord> >();
            Newton<DMatrixd, DVectord> nlsolver(lsolver);

            Parameters params;
            params.verbose(false);
            params.linear_solver_verbose(false);
            nlsolver.set_parameters(params);

            DVectord x = values(10, 2);
            TestFunctionND_1<DMatrixd, DVectord> fun2(x.size().get(0));

            DVectord expected = values(x.size().get(0), 0.468919);
            nlsolver.solve(fun2, x);
            assert(approxeq(expected, x));
           // std::cout << "         End: petsc_newton_test_outInfo" << std::endl;
        }

        void petsc_sparse_newton_test()
        {
            using namespace utopia;

            // std::cout << "         Begin: petsc_sparse_newton_test" << std::endl;
            auto lsolver = std::make_shared< BiCGStab<DSMatrixd, DVectord> >();
            Newton<DSMatrixd, DVectord> nlsolver(lsolver);
            nlsolver.enable_differentiation_control(false);

            Parameters params;
            params.verbose(false);
            params.linear_solver_verbose(false);
            nlsolver.set_parameters(params);

            SimpleQuadraticFunction<DSMatrixd, DVectord> fun;

            DVectord x = values(10, 2.);
            DVectord expected = zeros(x.size());

            nlsolver.solve(fun, x);
            assert(approxeq(expected, x));
           // std::cout << "         End: petsc_sparse_newton_test" << std::endl;
        }

        void petsc_newton_test()
        {
            using namespace utopia;

            // std::cout << "         Begin: PETSC_NEWTON_test" << std::endl;
            auto lsolver = std::make_shared< BiCGStab<DMatrixd, DVectord> >();
            Newton<DMatrixd, DVectord> nlsolver(lsolver);
            nlsolver.enable_differentiation_control(false);

            Parameters params;
            params.atol(1e-15);
            params.rtol(1e-15);
            params.stol(1e-15);
            params.verbose(false);
            nlsolver.set_parameters(params);

            SimpleQuadraticFunction<DMatrixd, DVectord> fun;

            DVectord x = values(_n, 2.);
            DVectord expected = zeros(x.size());

            nlsolver.solve(fun, x);
            assert(approxeq(expected, x));

            x = values(_n, 2.0);
            TestFunctionND_1<DMatrixd, DVectord> fun2(x.size().get(0));

            expected = values(x.size().get(0), 0.468919);
            nlsolver.solve(fun2, x);
            assert(approxeq(expected, x));

           // std::cout << "         End: PETSC_NEWTON_test" << std::endl;

        }

        void petsc_newton_rosenbrock_test()
        {
            using namespace utopia;

            auto lsolver = std::make_shared< BiCGStab<DMatrixd, DVectord> >();
            Newton<DMatrixd, DVectord> nlsolver(lsolver);
            nlsolver.enable_differentiation_control(false);

            Parameters params;
            params.rtol(1e-15);
            params.stol(1e-15);
            params.atol(1e-15);
            params.verbose(false);
            nlsolver.set_parameters(params);

            DVectord expected_rosenbrock;
            DVectord x0;

            if (mpi_world_size() <= 2) {
                RosenbrockGeneric<DMatrixd, DVectord> r_generic_2d;
                expected_rosenbrock = values(2, 1.0);
                x0 = values(2, 2.0);
                nlsolver.solve(r_generic_2d, x0);
                assert(approxeq(expected_rosenbrock, x0));
            }

            if (mpi_world_size() <= 3) {
                RosenbrockGeneric<DMatrixd, DVectord> r_generic_3d;
                expected_rosenbrock = values(3, 1.0);
                x0 = values(3, -2.0);
                nlsolver.solve(r_generic_3d, x0);
                assert(approxeq(expected_rosenbrock, x0));
            }

            if (mpi_world_size() <= 6) {
                RosenbrockGeneric<DMatrixd, DVectord> r_generic_6d;
                expected_rosenbrock = values(6, 1.0);
                x0 = values(6, 2.0);
                nlsolver.solve(r_generic_6d, x0);
                assert(approxeq(expected_rosenbrock, x0));
            }
        }

        void petsc_sparse_semismooth_newton_test()
        {
            using namespace utopia;

            // std::cout << "         Begin: petsc_sparse_semismooth_newton_test" << std::endl;
            auto lsolver = std::make_shared<BiCGStab<DSMatrixd, DVectord>>();
            DSMatrixd A, B;
            DVectord b, g;

            SemismoothNewton<DSMatrixd, DVectord> nlsolver(lsolver);
            nlsolver.enable_differentiation_control(false);

            Parameters params;
            params.verbose(false);
            nlsolver.set_parameters(params);

            // initial guess
            DVectord x_0 = values(_n, 60);
            {
                Write<DVectord> w(x_0);
                Range x_range = range(x_0);
                if (x_range.begin() == 0) x_0.set(0, 0);
                if (x_range.end() == _n) x_0.set(_n - 1, 0);
            }

            ExampleTestCase2<DSMatrixd, DVectord> example;
            example.getOperators(_n, A, b, g);
            nlsolver.solve(x_0, A, b, g);

           // std::cout << "         End: petsc_sparse_semismooth_newton_test" << std::endl;

        }

        void petsc_sparse_nonlinear_semismooth_newton_test()
        {
            using namespace utopia;

            // std::cout << "         Begin: petsc_sparse_nonlinear_semismooth_newton_test" << std::endl;
            auto lsolver = std::make_shared< ConjugateGradient<DSMatrixd, DVectord> >();

            NonlinSemismoothNewton<DSMatrixd, DVectord> nlsolver(lsolver);
            nlsolver.enable_differentiation_control(false);
            Parameters params;
            params.verbose(false);
            nlsolver.set_parameters(params);

            DSMatrixd A, B;
            DVectord upbo;

            ExampleTestCase<DSMatrixd, DVectord> example;
            example.getOperators(_n, A, B, upbo);

            DVectord rhs = values(_n, 60);
            {
                Write<DVectord> w(rhs);
                Range rhs_range = range(rhs);
                if (rhs_range.begin() == 0) rhs.set(0, 0);
                if (rhs_range.end() == _n) rhs.set(_n - 1, 0);

            }

            const QuadraticFunctionConstrained<DSMatrixd, DVectord> funn(rhs, A, B, upbo);
            nlsolver.solve(funn, rhs);
           // std::cout << "         End: petsc_sparse_nonlinear_semismooth_newton_test" << std::endl;

        }

        void petsc_direct_solver_newton_test()
        {
            using namespace utopia;

            // std::cout << "         Begin: petsc_direct_solver_newton_test" << std::endl;
            auto lsolver = std::make_shared< Factorization<DSMatrixd, DVectord> >();

#ifdef PETSC_HAVE_MUMPS
            lsolver->set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
#endif //PETSC_HAVE_MUMPS

            Newton<DSMatrixd, DVectord> nlsolver(lsolver);
            nlsolver.enable_differentiation_control(false);

            Parameters params;
            params.verbose(false);
            params.linear_solver_verbose(false);
            nlsolver.set_parameters(params);

            SimpleQuadraticFunction<DSMatrixd, DVectord> fun;

            DVectord x = values(_n, 2.);
            DVectord expected = zeros(x.size());

            nlsolver.solve(fun, x);
            assert(approxeq(expected, x));


            auto lCG = std::make_shared< ConjugateGradient<DSMatrixd, DVectord> >();
            nlsolver.set_linear_solver(lCG);
            x = values(_n, 2.);
            nlsolver.solve(fun, x);
            assert(approxeq(expected, x));
           // std::cout << "         End: petsc_direct_solver_newton_test" << std::endl;
        }

        void MG_test()
        {

            using namespace utopia;
            std::cout << "         Begin: MG_test" << std::endl;

            // reading data from outside
            DVectord rhs;
            DSMatrixd A, I_1, I_2, I_3;

            const std::string data_path = Utopia::Instance().get("data_path");

            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            // read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);


            std::vector <DSMatrixd> interpolation_operators;

            // from coarse to fine
            // interpolation_operators.push_back(std::move(I_1));
            interpolation_operators.push_back(std::move(I_2));
            interpolation_operators.push_back(std::move(I_3));

            //  init
            auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
#ifdef PETSC_HAVE_MUMPS
            direct_solver->set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
#endif //PETSC_HAVE_MUMPS

            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
            // auto smoother = std::make_shared<PointJacobi<DSMatrixd, DVectord>>();
            Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);


            multigrid.init_transfer(std::move(interpolation_operators));
            multigrid.galerkin_assembly(A);

            DVectord x_0 = zeros(A.size().get(0));

            Parameters params;
            params.verbose(false);
            params.linear_solver_verbose(false);
            multigrid.set_parameters(params);

            multigrid.solve(rhs, x_0);
           std::cout << "         End: MG_test" << std::endl;
        }


        void CG_MG_test()
        {
            // std::cout << "         Begin: CG_MG_test" << std::endl;

            //! [MG solve example]
            using namespace utopia;

            const bool verbose = false;

            DVectord rhs;
            DSMatrixd A, I_1, I_2, I_3;

            //reading data from disk
            const std::string data_path = Utopia::Instance().get("data_path");
            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            std::vector <DSMatrixd> interpolation_operators;

             //interpolation operators from coarse to fine
            interpolation_operators.push_back(std::move(I_1));
            interpolation_operators.push_back(std::move(I_2));
            interpolation_operators.push_back(std::move(I_3));

            //choose solver for coarse level solution
            auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();

#ifdef PETSC_HAVE_MUMPS
            direct_solver->set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
#endif //PETSC_HAVE_MUMPS

            //choose smoother
            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
            // auto smoother = std::make_shared<PointJacobi<DSMatrixd, DVectord>>();
            Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);
            multigrid.init_transfer(std::move(interpolation_operators));
            multigrid.max_it(1);
            multigrid.mg_type(2);


            ConjugateGradient<DSMatrixd, DVectord> cg;
            cg.verbose(verbose);

            DVectord x_0 = zeros(A.size().get(0));

            //CG with diagonal preconditioner
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<DSMatrixd, DVectord> >());
            cg.solve(A, rhs, x_0);

            //CG with multigrid preconditioner
            x_0 = zeros(A.size().get(0));
            cg.set_preconditioner(make_ref(multigrid));
            cg.solve(A, rhs, x_0);

            //Multigrid only
            // x_0 = zeros(A.size().get(0));
            // multigrid.max_it(12);
            // multigrid.verbose(verbose);
            // multigrid.solve(rhs, x_0);

            //! [MG solve example]


            // std::cout << "         End: CG_MG_test" << std::endl;
        }

        void PETSC_CG_MG_test()
        {
            using namespace utopia;

            const bool verbose = false;
            DVectord rhs;
            DSMatrixd A, I_1, I_2, I_3;

            const std::string data_path = Utopia::Instance().get("data_path");

            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            if(verbose) {
                disp("n unknowns:");
                disp(size(rhs));
            }

            std::vector <DSMatrixd> interpolation_operators;

             // from coarse to fine
            interpolation_operators.push_back(std::move(I_1));
            interpolation_operators.push_back(std::move(I_2));
            interpolation_operators.push_back(std::move(I_3));

            //  init
            auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
#ifdef PETSC_HAVE_SUPERLU_DIST
            direct_solver->set_type(SUPERLU_DIST_TAG, LU_DECOMPOSITION_TAG);
#else
            if(mpi_world_size() > 1) {
                if(mpi_world_rank() == 0) {
                    std::cerr << "[Error] Direct solver does not work in parallel compile with SuperLU" << std::endl;
                }
                return;
            }            
#endif //PETSC_HAVE_SUPERLU_DIST

            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
            Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);
            multigrid.init_transfer(std::move(interpolation_operators));
            multigrid.galerkin_assembly(A);

            multigrid.max_it(1);
            multigrid.mg_type(2);

            DVectord x_0;

            //! [KSPSolver solve example1]
            KSPSolver<DSMatrixd, DVectord> utopia_ksp;
//            GMRES<DSMatrixd, DVectord> utopia_ksp;
            auto precond = std::make_shared< InvDiagPreconditioner<DSMatrixd, DVectord> >();
            // utopia_ksp.set_preconditioner(precond);
            utopia_ksp.verbose(verbose);
            utopia_ksp.solve(A, rhs, x_0);
            //! [KSPSolver solve example1]

            x_0 = zeros(A.size().get(0));
            multigrid.verbose(false);
            multigrid.max_it(1);
            multigrid.mg_type(1);
            multigrid.pre_smoothing_steps(1);
            multigrid.post_smoothing_steps(1);
            utopia_ksp.set_preconditioner(make_ref(multigrid));
            utopia_ksp.solve(A, rhs, x_0);
        }


        void petsc_newton_PETScCG_utopia_test()
        {
            using namespace utopia;
            using namespace std;

            // std::cout << "         Begin: petsc_newton_PETScCG_utopia_test" << std::endl;

            //CG with diagonal preconditioner
            auto linear_solver  = make_shared< KSPSolver<DMatrixd, DVectord> >();
            auto preconditioner = make_shared< InvDiagPreconditioner<DMatrixd, DVectord> >();
            linear_solver->set_preconditioner(preconditioner);

            //Newton solver with cg linear solver
            Newton<DMatrixd, DVectord> newton_solver(linear_solver);
            newton_solver.verbose(false);

            const int n = 10;
            DVectord actual   = values(n, 2.);
            DVectord expected = values(n, 0.468919);

            TestFunctionND_1<DMatrixd, DVectord> fun(n);

            newton_solver.solve(fun, actual);
            assert(approxeq(expected, actual));
            // // std::cout << "         End: petsc_newton_PETScCG_utopia_test" << std::endl;
        }


        void petsc_newton_inexact_newton_with_KSP_test()
        {
            using namespace utopia;
            using namespace std;

            const bool verbose = false;

            auto linear_solver  = make_shared< KSPSolver<DMatrixd, DVectord> >();
            auto preconditioner = make_shared< InvDiagPreconditioner<DMatrixd, DVectord> >();
            linear_solver->set_preconditioner(preconditioner);


            InexactNewton<DMatrixd, DVectord> newton_solver(linear_solver);
            newton_solver.verbose(verbose);

            const int n = 10;
            DVectord actual   = values(n, 2.);
            DVectord expected = values(n, 0.468919);

            TestFunctionND_1<DMatrixd, DVectord> fun(n);

            newton_solver.solve(fun, actual);
            assert(approxeq(expected, actual));
            // // std::cout << "         End: petsc_newton_PETScCG_utopia_test" << std::endl;
        }





            PETScSolverTest()
        : _n(10) { }

    private:
        int _n;

    };
#endif //WITH_PETSC

    void runSolversTest()
    {
        std::cout << "Begin: SolverTest" << std::endl;
#ifdef WITH_PETSC
        SolverTest<DMatrixd, DVectord, PetscScalar>().run();
        PETScSolverTest().run();
#endif

#ifdef WITH_BLAS
        SolverTest<Matrixd, Vectord, double>().run();
#endif //WITH_BLAS

        std::cout << "End:   SolverTest" << std::endl;
    }
}

#include "utopia_Testing.hpp"
#include "utopia.hpp"

#ifdef WITH_PETSC

#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include "utopia_QuadraticFunction.hpp"
#include "utopia_petsc_TaoSolver.hpp"

namespace utopia {
    void petsc_tao_solve_simple()
    {
        TestFunctionND_1<PetscMatrix, PetscVector> fun(10);
        TaoSolver<PetscMatrix, PetscVector> tao(std::make_shared<Factorization<PetscMatrix, PetscVector>>());
        PetscVector x(layout(PetscCommunicator::get_default(), PetscTraits::decide(), 10), 0.0);
        tao.set_type("blmvm");
        tao.solve(fun, x);

        PetscVector expected(layout(x), 0.468919);
        utopia_test_assert(approxeq(x, expected));
    }

    void petsc_tao_solve_vi()
    {
        const SizeType n = 100;

        PetscMatrix m;
        PetscVector rhs, upper_bound;
        Poisson1D<PetscMatrix, PetscVector> ex2(n, 2);
        PetscVector x = ex2.initial_guess();
        ex2.hessian(x, m);
        ex2.get_rhs(rhs);
        upper_bound = ex2.upper_bound();

        auto lsolver = std::make_shared<ConjugateGradient<PetscMatrix, PetscVector>>();

        const double scale_factor = 1e-10;
        rhs *= scale_factor;
        upper_bound *= scale_factor;

        auto box = make_upper_bound_constraints(make_ref(upper_bound));

        QuadraticFunction<PetscMatrix, PetscVector> fun(make_ref(m), make_ref(rhs));
        TaoSolver<PetscMatrix, PetscVector> tao(lsolver);
        // tao.set_ksp_types("bcgs", "jacobi", " ");
        tao.set_box_constraints(box);
        // tao.set_type("tron");
        // tao.set_type("gpcg");
        tao.solve(fun, x);

        x *= 1./scale_factor;

        PetscVector xssn(layout(PetscCommunicator::get_default(), PetscTraits::decide(), n), 0.0);
        SemismoothNewton<PetscMatrix, PetscVector, HOMEMADE> ssnewton(std::make_shared<Factorization<PetscMatrix, PetscVector>>());
        ssnewton.set_box_constraints(box);
        ssnewton.stol(1e-18);
        ssnewton.atol(1e-18);
        ssnewton.rtol(1e-18);
        ssnewton.solve(m, rhs, xssn);

        xssn *= 1./scale_factor;

        double n_diff = norm2(xssn - x);
        utopia_test_assert(n_diff < 1e-10);


    }

    void petsc_tao_solve_mg()
    {
        PetscVector rhs;
        PetscMatrix A, I_1, I_2, I_3;

        const std::string data_path = Utopia::instance().get("data_path");

        read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
        read(data_path + "/laplace/matrices_for_petsc/f_A", A);
        read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
        read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

        std::vector<std::shared_ptr<PetscMatrix>> interpolation_operators;
        interpolation_operators.push_back(make_ref(I_2));
        interpolation_operators.push_back(make_ref(I_3));

        auto smoother      = std::make_shared<GaussSeidel<PetscMatrix, PetscVector>>();
        auto linear_solver = std::make_shared<ConjugateGradient<PetscMatrix, PetscVector>>();
        Multigrid<PetscMatrix, PetscVector> multigrid(smoother, linear_solver);
        multigrid.set_transfer_operators(std::move(interpolation_operators));
        PetscVector x(row_layout(A), 0.0);
        // PetscVector upper_bound = values(A.size().get(0), 0.003);
        // auto box = make_upper_bound_constraints(make_ref(upper_bound));

        QuadraticFunction<PetscMatrix, PetscVector> fun(make_ref(A), make_ref(rhs));
        TaoSolver<PetscMatrix, PetscVector> tao(make_ref(multigrid));

        // multigrid.verbose(true);
        multigrid.max_it(20);
        multigrid.atol(1e-15);
        multigrid.stol(1e-15);
        multigrid.rtol(1e-15);
        //constraints do not work with mg because system contains lagr mult
        // tao.set_box_constraints(box);
        tao.solve(fun, x);
    }

    void petsc_tao_tr_bound()
    {
        const SizeType n = 100;

        PetscMatrix m;
        PetscVector rhs, upper_bound;
        Poisson1D<PetscMatrix, PetscVector> ex2(n, 2);
        PetscVector x = ex2.initial_guess();
        ex2.hessian(x, m);
        ex2.get_rhs(rhs);
        upper_bound = ex2.upper_bound();

        const double scale_factor = 10e-10;
        rhs *= scale_factor;
        upper_bound *= scale_factor;

        auto box = make_upper_bound_constraints(make_ref(upper_bound));
        QuadraticFunction<PetscMatrix, PetscVector> fun(make_ref(m), make_ref(rhs));

        auto lsolver = std::make_shared<LUDecomposition<PetscMatrix, PetscVector> >();
        // auto lsolver = std::make_shared<BiCGStab<PetscMatrix, PetscVector> >();
        auto qp_solver = std::make_shared<TaoQPSolver<PetscMatrix, PetscVector> >(lsolver);

        // lsolver->atol(1e-16);

        qp_solver->atol(1e-15);
        qp_solver->stol(1e-15);
        qp_solver->max_it(10000);

        TrustRegionVariableBound<PetscMatrix, PetscVector>  tr_solver(qp_solver);
        tr_solver.set_box_constraints(box);
        tr_solver.verbose(false);
        tr_solver.atol(1e-15);
        tr_solver.stol(1e-15);
        tr_solver.rtol(1e-15);
        tr_solver.solve(fun, x);

        x *= 1./scale_factor;

        PetscVector xssn(layout(rhs), 0.0);
        SemismoothNewton<PetscMatrix, PetscVector, HOMEMADE> ssnewton(std::make_shared<Factorization<PetscMatrix, PetscVector>>());
        ssnewton.set_box_constraints(box);
        ssnewton.stol(1e-17);
        ssnewton.atol(1e-17);
        ssnewton.rtol(1e-17);
        ssnewton.solve(m, rhs, xssn);
        xssn *= 1./scale_factor;

        double n_diff = norm2(xssn - x);
        utopia_test_assert(n_diff < 1e-7);

    }


    static void tao()
    {
        //does not work yet missing ksp for dense matrix
        // UTOPIA_RUN_TEST(petsc_tao_solve_simple);
        // UTOPIA_RUN_TEST(petsc_tao_solve_vi);
        UTOPIA_RUN_TEST(petsc_tao_solve_mg);

//FIXME
#ifdef PETSC_HAVE_MUMPS
        UTOPIA_RUN_TEST(petsc_tao_tr_bound);
#endif //PETSC_HAVE_MUMPS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(tao);
}


#endif //WITH_PETSC

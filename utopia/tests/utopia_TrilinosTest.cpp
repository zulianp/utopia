
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS

#include "utopia_TrilinosTest.hpp"
#include "utopia.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_solvers.hpp"

#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include "test_problems/utopia_MultiLevelTestProblem.hpp"
#include <algorithm>

#ifdef WITH_PETSC
#include "utopia_petsc_trilinos.hpp"
#endif

#include "utopia_Structure.hpp"
#include "utopia_Eval_Structure.hpp"

#include "test_problems/utopia_BratuMultilevelTestProblem.hpp"
#include "test_problems/utopia_TestProblems.hpp"

#include "utopia_IPTransfer.hpp"

namespace utopia {

    template<class Matrix>
    static void build_rectangular_matrix(const SizeType &n, const SizeType &m, Matrix &mat)
    {
        mat  = local_sparse(n, m, 2);

        Write<TSMatrixd> w_(mat);
        auto r = row_range(mat);
        auto cols = size(mat).get(1);
        for(auto i = r.begin(); i < r.end(); ++i) {
            if(i >= cols) {
                break;
            }

            mat.set(i, i, 1.);

        }
    }

    template<class Matrix>
    static void build_rectangular_matrix_2(const SizeType &n, const SizeType &m, Matrix &mat)
    {
        mat  = local_sparse(n, m, 2);

        Write<TSMatrixd> w_(mat);
        auto r = row_range(mat);
        auto cols = size(mat).get(1);
        for(auto i = r.begin(); i < r.end(); ++i) {
            if(i >= cols) {
                break;
            }

            mat.set(i, i, 1.);

        }

        mat.set(0, m-1, 1.);
    }


    void trilinos_build()
    {
        auto n = 10;
        TVectord v = local_values(n, 1.);

        //FIXME replace this with an actual test
        // disp(v);

        TSMatrixd m = local_sparse(n, n, 3);
        assemble_laplacian_1D(m);

        //FIXME replace this with an actual test
        // disp(m);

        TVectord z = m * v;
        double nz = norm2(z);
        utopia_test_assert(approxeq(nz, 0.));
    }

    void trilinos_build_identity()
    {
        auto n = 10;
        TSMatrixd id = local_identity(n, n);
        TVectord v = local_values(n, 2.);
        double actual = norm1(id * v);

        utopia_test_assert(approxeq(size(v).get(0) * 2., actual));

        TSMatrixd id_t = transpose(id);
        actual = norm1(id * v);

        utopia_test_assert(approxeq(size(v).get(0) * 2., actual));
    }

    void trilinos_rect_matrix()
    {
        TSMatrixd P;
        // build_rectangular_matrix(5, 10, P);
        build_rectangular_matrix(10, 5, P);


        // auto rm = P.implementation().implementation().getRangeMap();
        // auto dm = P.implementation().implementation().getDomainMap();

        // auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

        // rm->describe(*out);
        // dm->describe(*out);

        utopia_test_assert(P.implementation().is_valid(true));
    }

    void trilinos_accessors()
    {
        auto n = 10;
        TVectord v = local_values(n, 10.);

        {
            Range r = range(v);
            Write<TVectord> w_v(v);

            //set first and last entries of each process, are to be 0
            v.set(r.begin(), 0.);
            v.add(r.end() - 1, -10.);
        }

        Size s = size(v);
        Size ls = local_size(v);

        utopia_test_assert(s.get(0) == n * mpi_world_size());
        utopia_test_assert(ls.get(0) == n);

        //FIXME replace this with an actual test
        // disp(v);
    }

    void trilinos_vec_axpy()
    {
        auto n = 10;
        TVectord y = local_values(n, 1.);
        TVectord x = local_values(n, 5.);
        auto alpha = 0.1;
        y += alpha * x;

        double val = norm1(y);
        utopia_test_assert(approxeq(val, n * mpi_world_size() * 1.5));
    }

    void trilinos_vec_scale()
    {
        auto n = 10;
        TVectord y = local_values(n, 1.);

        y *= 2.;

        double val = norm1(y);
        utopia_test_assert(approxeq(val, size(y).get(0) * 2.));


        TVectord y2 = y * 2.;

        val = norm1(y2);
        utopia_test_assert(approxeq(val, size(y2).get(0) * 4.));
    }

    void trilinos_mat_scale()
    {
        auto n = 10;
        TSMatrixd Y = local_sparse(n, n, 3);
        assemble_laplacian_1D(Y);

        Y *= 2.;

        TVectord x = local_values(n, 1.);
        double val = norm1(Y * x);
        utopia_test_assert(approxeq(val, 0.));


        TSMatrixd Y2 = Y * 2.;

        val = norm1(Y2 * x);
        utopia_test_assert(approxeq(val, 0.));
    }

    void trilinos_mat_axpy()
    {
        auto n = 10;
        TSMatrixd X = local_sparse(n, n, 3);
        assemble_laplacian_1D(X);

        TSMatrixd Y = X;

        auto alpha = 0.1;
        Y += alpha * X;


        TVectord v = local_values(n, 5.);

        double val = norm1(Y * v);
        utopia_test_assert(approxeq(val, 0.));
    }

    void trilinos_mv()
    {
        auto n = 10;
        TVectord x = local_values(n, 5.);
        TSMatrixd m = sparse(n, n, 3);
        assemble_laplacian_1D(m);

        TVectord y = m * x;

        const double val = norm2(y);
        utopia_test_assert(approxeq(val, 0.));
    }

    void trilinos_transpose()
    {
        auto rows = 5;
        auto cols = 11;
        TSMatrixd A = local_sparse(rows, cols, 2);

        {
            Write<TSMatrixd> w_A(A);
            Range r = row_range(A);

            for(auto i = r.begin(); i < r.end(); ++i) {
                A.set(i, 0, 1.);
                A.set(i, 1, 2.);
            }
        }

        auto s = size(A);

        TSMatrixd At = transpose(A);
        auto s_t = size(At);
        utopia_test_assert(s_t.get(0) == s.get(1));


        TSMatrixd id = local_identity(rows, cols);
        TVectord v = local_values(cols, 2.);

        TVectord actual = id * v;
        double norm_actual = norm1(actual);

        utopia_test_assert(approxeq(size(actual).get(0) * 2, norm_actual));

        //Does not work in parallel
        TSMatrixd id_t = transpose(id);
        TVectord v2 = local_values(rows, 2.);
        actual = id_t * v2;
        norm_actual = norm1(actual);
        double norm_expected = size(v2).get(0) * 2.;

        // disp(id);
        // std::cout << "-----------------------" << std::endl;
        // disp(id_t);
        // std::cout << "-----------------------" << std::endl;

        // std::cout << norm_expected << " == " << norm_actual << std::endl;
        utopia_test_assert(approxeq(norm_expected, norm_actual));
    }


    void trilinos_mm()
    {
        auto n = 10;
        TSMatrixd A = local_sparse(n, n, 3);
        assemble_laplacian_1D(A);
        TSMatrixd B = A;
        TSMatrixd C = transpose(A) * B;

        TVectord x = local_values(n, 5.);
        TVectord y = C * x;

        const double val = norm2(y);
        utopia_test_assert(approxeq(val, 0.));
    }

    void trilinos_m_tm()
    {
        auto n = 10;
        auto m = 3;
        TSMatrixd A = local_sparse(n, n, 3);
        assemble_laplacian_1D(A);
        TSMatrixd P;
        build_rectangular_matrix(n, m, P);
        TSMatrixd B = A * P;
        TSMatrixd P_t = transpose(P);
        TSMatrixd C_1 = P_t * A;
        TSMatrixd C_2 = transpose(P) * A;


        //FIXME write test here
    }

    void trilinos_diag()
    {
        auto n = 10;
        TSMatrixd A = local_sparse(n, n, 3);
        assemble_laplacian_1D(A);
        TVectord d = diag(A);

        const double val = norm1(d);
        utopia_test_assert(approxeq(val, size(d).get(0)*2.-2.));

        TSMatrixd D = diag(d);
        TVectord x  = local_values(n, 1.);
        utopia_test_assert(approxeq(d, D*x));
    }

    void trilinos_ptap()
    {
        auto n = 10;
        auto m = 3;
        TSMatrixd A = local_sparse(n, n, 3);
        assemble_laplacian_1D(A);

        TSMatrixd P;
        build_rectangular_matrix(n, m, P);

        //For the moment this is computing (transpose(P) * A) * P
        TSMatrixd R   = utopia::ptap(A, P);
        //same thing
        TSMatrixd R_2 = transpose(P) * A * P;

        // disp(A);
        // disp(P);
        // disp(R);

        //FIXME write test here
    }

    void trilinos_cg()
    {
         MultiLevelTestProblem<TSMatrixd, TVectord> ml_problem(100, 2);
         TVectord x = zeros(size(*ml_problem.rhs));
         (*ml_problem.rhs) *= 0.0001;

         ConjugateGradient<TSMatrixd, TVectord> cg;
         cg.rtol(1e-6);
         cg.atol(1e-6);
         cg.max_it(500);
         // cg.verbose(true);
         cg.update(ml_problem.matrix);
         cg.apply(*ml_problem.rhs, x);

         utopia_test_assert(approxeq(*ml_problem.rhs, *ml_problem.matrix * x, 1e-5));
    }

    template<class Matrix, class Vector>
    void test_mg()
    {
        using TransferT       = utopia::Transfer<Matrix, Vector>;
        using IPTransferT     = utopia::IPTransfer<Matrix, Vector>;
        using MatrixTransferT = utopia::MatrixTransfer<Matrix, Vector>;

        const static bool verbose   = false;
        const static bool use_masks = true;

        MultiLevelTestProblem<Matrix, Vector> ml_problem(5, 2, !use_masks);

        auto smoother      = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
        auto coarse_solver = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();

        Multigrid<Matrix, Vector> multigrid(
            smoother,
            coarse_solver
        );

        multigrid.max_it(1);
        multigrid.atol(1e-13);
        multigrid.stol(1e-13);
        multigrid.rtol(1e-9);
        multigrid.pre_smoothing_steps(3);
        multigrid.post_smoothing_steps(3);
        multigrid.set_fix_semidefinite_operators(true);
        multigrid.must_generate_masks(use_masks);;
        multigrid.verbose(verbose);

        std::vector<std::shared_ptr<TransferT>> transfers;

        for(auto &interp_ptr : ml_problem.interpolators) {
            if(use_masks) {
                //compute transpose explicitly for restriction
                transfers.push_back( std::make_shared<MatrixTransferT>(interp_ptr) );
            } else {
                //apply transpose for restriction
                transfers.push_back( std::make_shared<IPTransferT>(interp_ptr) );
            }

            // utopia_test_assert(interp_ptr->implementation().is_valid(true));
        }

        multigrid.set_transfers(transfers);

        Vector x = zeros(size(*ml_problem.rhs));
        multigrid.update(ml_problem.matrix);

        if(verbose) {
            multigrid.describe();
        }

        multigrid.apply(*ml_problem.rhs, x);

        double diff0 = norm2(*ml_problem.matrix * x);
        double diff  = norm2(*ml_problem.rhs - *ml_problem.matrix * x);
        double rel_diff = diff/diff0;

        // utopia_test_assert(rel_diff < 1e-8);
    }

    void trilinos_mg_1D()
    {
        // if(mpi_world_size() > 1) return;
      //petsc version
      test_mg<DSMatrixd, DVectord>();
      
      std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
      //trilinos version
      test_mg<TSMatrixd, TVectord>();
    }


    void trilinos_mg()
    {
        // if(mpi_world_size() > 1) return;

        bool ok = true;

        TVectord rhs;
        TSMatrixd A, I;

        Multigrid<TSMatrixd, TVectord> multigrid(
            std::make_shared<ConjugateGradient<TSMatrixd, TVectord>>(),
            std::make_shared<ConjugateGradient<TSMatrixd, TVectord>>()
        );


#ifdef WITH_PETSC
        //FIXME needs trilinos formats but for the moment lets use petsc's
        {
            DSMatrixd petsc_A, petsc_I;
            DVectord petsc_rhs;

            const std::string folder =  Utopia::instance().get("data_path") + "/mg";

            ok = read(folder + "/rhs.bin", petsc_rhs); utopia_test_assert(ok);
            ok = read(folder + "/A.bin", petsc_A);     utopia_test_assert(ok);
            ok = read(folder + "/I.bin", petsc_I);     utopia_test_assert(ok);

            backend_convert_sparse(petsc_I, I);
            backend_convert_sparse(petsc_A, A);
            backend_convert(petsc_rhs, rhs);
        }

        // write("A.mm", A);
        // write("I.mm", I);

        std::vector<std::shared_ptr<TSMatrixd>> interpolation_operators;
        interpolation_operators.push_back(make_ref(I));

        multigrid.set_transfer_operators(std::move(interpolation_operators));
        multigrid.max_it(20);
        multigrid.atol(1e-15);
        multigrid.stol(1e-15);
        multigrid.rtol(1e-15);
        multigrid.verbose(true);
        multigrid.set_fix_semidefinite_operators(true);
        multigrid.must_generate_masks(true);
        TVectord x = local_zeros(local_size(rhs));

        try {
            multigrid.update(make_ref(A));
            ok = multigrid.apply(rhs, x); utopia_test_assert(ok);
        } catch(const std::exception &ex) {
            std::cout << ex.what() << std::endl;
            utopia_test_assert(false);
        }

        std::cout << std::flush;

        double diff = norm2(rhs - A * x);
        utopia_test_assert(approxeq(diff, 0., 1e-6));

#endif //WITH_PETSC

    }

    void trilinos_row_view()
    {
        TSMatrixd A = local_sparse(4, 4, 3);
        assemble_laplacian_1D(A);

        auto rr = row_range(A);

        for(auto i = rr.begin(); i != rr.end(); ++i) {
            RowView<TSMatrixd> row(A, i);
            utopia_test_assert(row.n_values() >= 2);
            auto col = row.col(0);

            utopia_test_assert(col == i || col == i - 1 || col == i  + 1);
        }
    }

    void trilinos_row_view_and_loops()
    {
        int n = 10;
        int m = 3;


        TSMatrixd P;
        build_rectangular_matrix(n, m, P);

        auto rr = row_range(P);

        int nnz = std::min(n, int(size(P).get(1)));
        if(rr.begin() >= nnz)
        {
            nnz = 0;
        }

        SizeType count = 0;
        each_read(P, [&count](const SizeType i, const SizeType j, const double val) {
            utopia_test_assert(val == 1.);
            ++count;
        });

        utopia_test_assert(nnz == count);

        TSMatrixd P_t = transpose(P);

        each_read(P_t, [&count](const SizeType i, const SizeType j, const double val) {
            utopia_test_assert(val == 1.);
            --count;
        });

        if(mpi_world_size() == 1) {
            utopia_test_assert(count == 0);
        }

        // disp(P);
        // disp(P_t);
    }



    void trilinos_each_read_transpose()
    {
        MultiLevelTestProblem<TSMatrixd, TVectord> ml_problem(5, 2, false);


        TSMatrixd R = transpose(*ml_problem.interpolators[0]);
        TSMatrixd R_copy = R;
        R_copy *= 0.;

        TVectord v  = local_values(local_size(R).get(1), 10.);
        TVectord Rv = local_zeros(local_size(R).get(0));

        Rv = R * v;

        double nrv = norm2(Rv);
        utopia_test_assert(nrv > 10.);
    }


    void trilinos_read()
    {
        TSMatrixd m;
        auto path = Utopia::instance().get("data_path") + "/matrixmarket/gre_343_343_crg.mm";
        bool ok = read(path, m);
        utopia_test_assert(ok);
    }


#ifdef WITH_PETSC
    void trilinos_petsc_interop()
    {
        KSPSolver<TSMatrixd, TVectord> solver;

        MultiLevelTestProblem<TSMatrixd, TVectord> ml_problem(10, 2);
        TVectord x = zeros(size(*ml_problem.rhs));
        (*ml_problem.rhs) *= 0.0001;


        solver.solve(*ml_problem.matrix, *ml_problem.rhs, x);

        // disp(*ml_problem.matrix);

        DSMatrixd p_mat;
        backend_convert_sparse(*ml_problem.matrix, p_mat);

        // disp(p_mat);
        double diff = norm2(*ml_problem.rhs - *ml_problem.matrix * x);
        utopia_test_assert(approxeq(diff, 0., 1e-8));
    }
#endif //WITH_PETSC

    void trilinos_structure()
    {
        auto n = 10;
        TSMatrixd A = local_sparse(n, n, 3);
        assemble_laplacian_1D(A);

        auto expr = structure(A);

        TSMatrixd B(expr);
    }

    void trilinos_rmtr()
    {
        BratuMultilevelTestProblem<TSMatrixd, TVectord> problem;
        problem.verbose = true;

        TVectord x = values(problem.n_dofs[problem.n_levels -1 ], 0.0);

        std::vector<std::shared_ptr<ExtendedFunction<TSMatrixd, TVectord> > >  level_functions(problem.n_levels);


        for(auto l=0; l < problem.n_levels; l++)
        {
            Bratu1D<TSMatrixd, TVectord> fun(problem.n_dofs[l]);
            level_functions[l] = std::make_shared<Bratu1D<TSMatrixd, TVectord> >(fun);

            // making sure that fine level IG is feasible
            if(l+1 == problem.n_levels)
                fun.apply_bc_to_initial_guess(x);
        }

        auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<TSMatrixd, TVectord> >();
        tr_strategy_coarse->atol(1e-12);
        tr_strategy_coarse->rtol(1e-12);

        auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<TSMatrixd, TVectord> >();
        tr_strategy_fine->atol(1e-12);
        tr_strategy_fine->rtol(1e-12);

        // auto rmtr = std::make_shared<RMTR<TSMatrixd, TVectord, SECOND_ORDER>  >(tr_strategy_coarse, tr_strategy_fine);
        auto rmtr = std::make_shared<RMTR<TSMatrixd, TVectord, GALERKIN>  >(tr_strategy_coarse, tr_strategy_fine);
        rmtr->set_transfer_operators(problem.prolongations, problem.restrictions);

        rmtr->max_it(1000);
        rmtr->max_coarse_it(1);
        rmtr->max_smoothing_it(1);
        rmtr->delta0(1);
        rmtr->atol(1e-6);
        rmtr->rtol(1e-10);
        rmtr->set_grad_smoothess_termination(0.000001);
        rmtr->set_eps_grad_termination(1e-7);

        rmtr->verbose(problem.verbose);
        // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
        rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
        rmtr->set_functions(level_functions);


        rmtr->solve(x);
    }

    void run_trilinos_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("TrilinosTest");
        UTOPIA_RUN_TEST(trilinos_structure);
        UTOPIA_RUN_TEST(trilinos_build);
        UTOPIA_RUN_TEST(trilinos_build_identity);
        UTOPIA_RUN_TEST(trilinos_accessors);
        UTOPIA_RUN_TEST(trilinos_vec_scale);
        UTOPIA_RUN_TEST(trilinos_mat_scale);
        UTOPIA_RUN_TEST(trilinos_vec_axpy);
        UTOPIA_RUN_TEST(trilinos_mat_axpy);
        UTOPIA_RUN_TEST(trilinos_mv);
        UTOPIA_RUN_TEST(trilinos_mm);
        UTOPIA_RUN_TEST(trilinos_m_tm);
        UTOPIA_RUN_TEST(trilinos_diag);
        UTOPIA_RUN_TEST(trilinos_read);
        UTOPIA_RUN_TEST(trilinos_ptap);
        UTOPIA_RUN_TEST(trilinos_cg);
        UTOPIA_RUN_TEST(trilinos_rect_matrix);


        //tests that fail in parallel
        UTOPIA_RUN_TEST(trilinos_mg_1D);
        UTOPIA_RUN_TEST(trilinos_row_view);
        UTOPIA_RUN_TEST(trilinos_row_view_and_loops);
        UTOPIA_RUN_TEST(trilinos_transpose);
        UTOPIA_RUN_TEST(trilinos_each_read_transpose);

#ifdef WITH_PETSC
        UTOPIA_RUN_TEST(trilinos_petsc_interop);
#endif //WITH_PETSC

        //tests that always fail
        UTOPIA_RUN_TEST(trilinos_rmtr);
        // UTOPIA_RUN_TEST(trilinos_mg);

        UTOPIA_UNIT_TEST_END("TrilinosTest");
    }
}

#else //WITH_TRILINOS

namespace utopia
{
    void run_trilinos_test() {}
}

#endif //WITH_TRILINOS

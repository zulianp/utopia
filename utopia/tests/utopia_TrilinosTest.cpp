#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS

#include "utopia_TrilinosTest.hpp"
#include "utopia.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_solvers.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include "test_problems/utopia_MultiLevelTestProblem.hpp"
#include <algorithm>

namespace utopia {   
    
    template<class TensorFrom, class TensorTo>
    void backend_convert_sparse(const Wrapper<TensorFrom, 2> &from, Wrapper<TensorTo, 2> &to)
    {
        auto ls = local_size(from);
        auto n_row_local = ls.get(0);
        std::vector<int> nnzxrow(n_row_local, 0);
        auto r = row_range(from);
        
        each_read(from, [&nnzxrow,&r](const SizeType i, const SizeType j, const double val) {
            ++nnzxrow[i - r.begin()];
        });
        
        
        auto nnz = *std::max_element(nnzxrow.begin(), nnzxrow.end());
        
        //FIXME use nnzxrow instead
        to = local_sparse(ls.get(0), ls.get(1), nnz);
        
        
        Write<Wrapper<TensorTo, 2>> w_t(to);
        each_read(from, [&to](const SizeType i, const SizeType j, const double val) {
            to.set(i, j, val);
        });
    }
    
    template<class TensorFrom, class TensorTo>
    void backend_convert(const Wrapper<TensorFrom, 1> &from, Wrapper<TensorTo, 1> &to)
    {
        auto ls = local_size(from).get(0);
        to = local_zeros(ls);
        
        Write< Wrapper<TensorTo, 1> > w_t(to);
        each_read(from, [&to](const SizeType i, const double val) {
            to.set(i, val);
        });
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
    
    template<class Matrix>
    static void build_rectangular_matrix(const SizeType &n, const SizeType &m, Matrix &mat)
    {
        mat  = local_sparse(n, m, 1);
        
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
    
    void trilinos_belos()
    {
        // comunicator
        // map
        
        TVectord x ;
        TVectord b ;
        TSMatrixd A;
        Parameters params;
        //params.set_param_file_name( "~/utopiaTrilinosFile.xml");
        
        BelosSolver<TSMatrixd, TVectord> solver(params);
        
        solver.solve(A, b, x);
        std::cout << "Number of Iterations " << solver.getNumIter() << std::endl;
        
        std::cout << "Achieved tolerance " << solver.achievedTol() << std::endl;
        
        /*   Parameters<TRILINOS> param;
         PrecondionedSolver<TSMatrixd, TVectord, TVectord, TRILINOS> prec;
         prec.set_preconditioner();
         
         linearSol.apply();
         ///////////////////
         
         Teuchos::RCP<const Teuchos::Comm<int> > Comm =
         Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
         int myPID = Comm->getRank();
         
         Teuchos::RCP<const map_type> Map =
         Teuchos::rcp(new map_type(Ndofs, indexBase, Comm));
         
         Teuchos::RCP<vec_type> LHS = Teuchos::rcp(new vec_type(Map, false));
         Teuchos::RCP<vec_type> RHS = Teuchos::rcp(new vec_type(Map, false));
         
         for (int i = 0; i < Ndofs; ++i)
         LHS->replaceLocalValue(i, _linSystem->solution(i));
         LHS->getData(0)[i];
         
         Teuchos::RCP<matrix_type> A = Teuchos::rcp(
         new matrix_type(Map, maxNumEntries, Tpetra::StaticProfile));
         
         A->insertGlobalValues(row, values.size(), values.data(), columns.data());
         A->fillComplete();
         
         Teuchos::RCP<Belos::LinearProblem<SC, mv_type, op_type> > linearProblem =
         Teuchos::rcp(new problem_type(A, LHS, RHS));
         
         linearProblem->setProblem();
         
         //list
         Teuchos::RCP<Teuchos::ParameterList> ParamList = Teuchos::getParametersFromXmlFile(param_file_name);
         //sublist
         auto& fasterPL = ParamList->sublist("FASTER", true);
         bool direct_solver = fasterPL.get("Direct Solver", false);
         bool direct_prec = fasterPL.get<bool>("Direct Preconditioner", false);
         std::string dir_prec_type = fasterPL.get("Ifpack2 Preconditioner", "prec_type_unset");
         std::string sol_type = fasterPL.get("Solver Type", "CG");
         //"factory"
         Teuchos::RCP<Amesos2::Solver<matrix_type, mv_type> > directSolver;
         Teuchos::RCP<solver_type> belosSolver;
         if (false==direct_solver)
         {//belos
         { Belos::SolverFactory<SC, mv_type, op_type> belosFactory;
         belosSolver = belosFactory.create(sol_type, Teuchos::sublist(ParamList, sol_type, false));
         }
         //preconditioner
         Teuchos::RCP<ifpack_prec_type> M_ifpack;
         Teuchos::RCP<muelu_prec_type> M_muelu;
         if (direct_prec) {
         M_ifpack = Ifpack2::Factory::create<matrix_type>(dir_prec_type, A);
         assert(!M_ifpack.is_null());
         M_ifpack->setParameters(ParamList->sublist(dir_prec_type, false));
         M_ifpack->initialize();
         M_ifpack->compute();
         linearProblem->setLeftPrec(M_ifpack);
         } else {
         // Multigrid Hierarchy
         M_muelu = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<op_type>)A,
         ParamList->sublist("MueLu", false));
         assert(!M_muelu.is_null());
         linearProblem->setRightPrec(M_muelu);
         }
         //solve
         belosSolver->setProblem(linearProblem);
         belosSolver->getCurrentParameters()->print();
         const Belos::ReturnType belosResult = belosSolver->solve();
         //print
         ret = belosResult;
         int numIterations = belosSolver->getNumIters();
         residualOut = belosSolver->achievedTol();
         if (myPID == 0) {
         std::cout << "number of iterations = " << numIterations << std::endl;
         std::cout << "||Residual|| = " << residualOut << std::endl;
         }
         }else//amesos
         {
         direct_solver = true;
         directSolver =
         Amesos2::create<matrix_type, mv_type>(sol_type, A, RHS, LHS);
         directSolver->setParameters(Teuchos::sublist(ParamList, sol_type, false));
         directSolver->symbolicFactorization().numericFactorization().solve();
         
         
         }*/
        ////////////////////
        
        
    }
    
    
    void trilinos_mg_1D()
    {
        if(mpi_world_size() > 1) return;
        
        const static bool verbose = true;
        
        MultiLevelTestProblem<TSMatrixd, TVectord> ml_problem(4, 2, false);
        // ml_problem.write_matlab("./");
        
        auto smoother = std::make_shared<ConjugateGradient<TSMatrixd, TVectord>>();
        auto coarse_solver = std::make_shared<ConjugateGradient<TSMatrixd, TVectord>>();
        Multigrid<TSMatrixd, TVectord> multigrid(
                                                 smoother,
                                                 coarse_solver
                                                 );
        
        // smoother->verbose(true);
        // coarse_solver->verbose(true);
        
        multigrid.set_transfer_operators(ml_problem.interpolators);
        multigrid.max_it(4);
        multigrid.atol(1e-12);
        multigrid.stol(1e-10);
        multigrid.rtol(1e-10);
        multigrid.pre_smoothing_steps(3);
        multigrid.post_smoothing_steps(3);
        multigrid.set_fix_semidefinite_operators(true);
        multigrid.verbose(verbose);
        
        TVectord x = zeros(size(*ml_problem.rhs));
        multigrid.update(ml_problem.matrix);
        
        // write("A0.txt", multigrid.level(0).A());
        // write("R0.txt", multigrid.transfer(0).R());
        
        if(verbose) {
            multigrid.describe();
        }
        
        multigrid.apply(*ml_problem.rhs, x);
        
        double diff = norm2(*ml_problem.rhs - *ml_problem.matrix * x);
        disp(diff);
        
        utopia_test_assert(approxeq(*ml_problem.rhs, *ml_problem.matrix * x, 1e-7));
    }
    
    void trilinos_mg()
    {
        // if(mpi_world_size() > 1) return;
        
        bool ok = true;
        
        TVectord rhs;
        TSMatrixd A, I_1, I_2, I_3;
        
        Multigrid<TSMatrixd, TVectord> multigrid(
                                                 std::make_shared<ConjugateGradient<TSMatrixd, TVectord>>(),
                                                 std::make_shared<ConjugateGradient<TSMatrixd, TVectord>>()
                                                 );
        
        
#ifdef WITH_PETSC        
        //FIXME needs trilinos formats but for the moment lets use petsc's
        {
            DSMatrixd petsc_A, petsc_I_1, petsc_I_2, petsc_I_3;
            DVectord petsc_rhs;
            
            const std::string folder =  Utopia::instance().get("data_path") + "/laplace/matrices_for_petsc";
            
            ok = read(folder + "/f_rhs", petsc_rhs); utopia_test_assert(ok);
            ok = read(folder + "/f_A", petsc_A);     utopia_test_assert(ok);
            ok = read(folder + "/I_1", petsc_I_1);   utopia_test_assert(ok);
            ok = read(folder + "/I_2", petsc_I_2);   utopia_test_assert(ok);
            ok = read(folder + "/I_3", petsc_I_3);   utopia_test_assert(ok);
            
            backend_convert_sparse(petsc_I_1, I_1);
            backend_convert_sparse(petsc_I_2, I_2);
            backend_convert_sparse(petsc_I_3, I_3);
            backend_convert_sparse(petsc_A, A);
            backend_convert(petsc_rhs, rhs);
        }
        
        std::vector<std::shared_ptr<TSMatrixd>> interpolation_operators;
        // interpolation_operators.push_back(make_ref(I_1));
        // interpolation_operators.push_back(make_ref(I_2));
        interpolation_operators.push_back(make_ref(I_3));
        
        multigrid.set_transfer_operators(std::move(interpolation_operators));
        multigrid.max_it(20);
        multigrid.atol(1e-15);
        multigrid.stol(1e-15);
        multigrid.rtol(1e-15);
        multigrid.verbose(true);
        multigrid.must_generate_masks(false);
        TVectord x = local_zeros(local_size(rhs));
        
        try {
            multigrid.update(make_ref(A));
            ok = multigrid.apply(rhs, x); utopia_test_assert(ok);
        } catch(const std::exception &ex) {
            std::cout << ex.what() << std::endl;
            utopia_test_assert(false);
        }
        
        std::cout << std::flush;
        utopia_test_assert(approxeq(rhs, A * x, 1e-6));
#endif //WITH_PETSC
        
    }
    
    void row_view_and_loops()
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
    
    void trilinos_read()
    {
        TSMatrixd m;
        auto path = Utopia::instance().get("data_path") + "/matrixmarket/gre_343_343_crg.mm";
        bool ok = read(path, m);
        utopia_test_assert(ok);
    }
    
    void run_trilinos_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("TrilinosTest");
        UTOPIA_RUN_TEST(trilinos_build);
        UTOPIA_RUN_TEST(trilinos_build_identity);
        UTOPIA_RUN_TEST(trilinos_accessors);
        UTOPIA_RUN_TEST(trilinos_vec_scale);
        UTOPIA_RUN_TEST(trilinos_mat_scale);
        UTOPIA_RUN_TEST(trilinos_vec_axpy);
        UTOPIA_RUN_TEST(trilinos_mat_axpy);
        UTOPIA_RUN_TEST(trilinos_transpose);
        UTOPIA_RUN_TEST(trilinos_mv);
        UTOPIA_RUN_TEST(trilinos_mm);
        UTOPIA_RUN_TEST(trilinos_m_tm);
        UTOPIA_RUN_TEST(trilinos_diag);
        UTOPIA_RUN_TEST(trilinos_read);
        UTOPIA_RUN_TEST(trilinos_ptap);
        UTOPIA_RUN_TEST(trilinos_cg);
        UTOPIA_RUN_TEST(trilinos_belos);
        
        //tests that fail in parallel
        UTOPIA_RUN_TEST(row_view_and_loops);
        
        //tests that always fail
        // UTOPIA_RUN_TEST(trilinos_mg_1D);
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

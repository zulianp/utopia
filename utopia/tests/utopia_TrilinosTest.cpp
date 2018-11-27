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

#include <cmath>

#include "utopia_trilinos_Each.hpp"

namespace utopia {
    
    template<class Matrix>
    static void build_rectangular_matrix(const SizeType &n, const SizeType &m, Matrix &mat)
    {
        mat  = local_sparse(n, m, 2);
        
        Write<Matrix> w_(mat);
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
        
        Write<Matrix> w_(mat);
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
    
    void trilinos_matrix_access()
    {
        auto n = 10;
        TSMatrixd Y = local_sparse(n, n, 3);
        assemble_laplacian_1D(Y);
        
        
        auto rr = row_range(Y);
        auto i  = rr.begin();
        
        {
            Read<TSMatrixd> r_(Y);
            if(i == 0 || i + 1 == size(Y).get(0)) {
                utopia_test_assert(approxeq(Y.get(i, i), 1.));
            } else {
                utopia_test_assert(approxeq(Y.get(i, i), 2.));
            }
        }
        
        {
            Write<TSMatrixd> w_(Y);
            Y.set(i, i, 4.);
        }
        
        {
            Read<TSMatrixd> r_(Y);
            utopia_test_assert(approxeq(Y.get(i, i), 4.));
        }
    }
    
    void trilinos_set()
    {
        auto n = 10;
        TVectord v = local_values(n, 10.);
        
        v.set(0.);
        
        double sum_v = sum(v);
        
        utopia_test_assert(approxeq(sum_v, 0.));
    }
    
    void trilinos_vec_minus()
    {
        auto n = 10;
        TVectord y = local_values(n, 1.);
        TVectord x = local_values(n, 5.);
        
        TVectord z;
        z = y - x;
        
        TVectord expected = local_values(n, -4.);
        
        double sum_z = double(sum(z));
        
        utopia_test_assert(approxeq(sum_z, size(z).get(0) * (-4.)));
        utopia_test_assert(approxeq(z, expected));
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
    
    void trilinos_residual()
    {
        auto n = 10;
        TVectord y   = local_values(n, 2.);
        TVectord x   = local_values(n, 1.);
        TSMatrixd Id = local_identity(n, n);
        
        TVectord z = x - Id * y;
        
        double val = norm1(z);
        utopia_test_assert(approxeq(val, size(y).get(0)));
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
        double tolerance = 30. * std::numeric_limits<double>::epsilon();
        //std::cout << "val " << val <<std::endl;
        utopia_test_assert(approxeq(val, 0., tolerance ));
        
        TSMatrixd Id = local_identity(n, n);
        Id += 2. * Id;
        
        v.set(1.);
        val = norm1(Id * v);
        utopia_test_assert( approxeq(val, size(v).get(0) * 3., 1e-14));
    }
    
    void trilinos_mv()
    {
        auto n = 10;
        TVectord x = local_values(n, 5.);
        TSMatrixd m = local_sparse(n, n, 3);
        assemble_laplacian_1D(m);
        
        TVectord y = m * x;
        
        const double val = norm2(y);
        utopia_test_assert(approxeq(val, 0.));
    }
    
    
    void trilinos_apply_transpose()
    {
        auto rows = 5;
        auto cols = 6;
        TSMatrixd A = local_sparse(rows, cols, 2);
        
        {
            Write<TSMatrixd> w_A(A);
            Range r = row_range(A);
            
            for(auto i = r.begin(); i < r.end(); ++i) {
                A.set(i, i,     1.);
                A.set(i, i + 1, 1.);
            }
        }
        
        TVectord v    = local_values(rows, 1.);
        TVectord At_v = transpose(A) * v;
        
        each_read(At_v, [](const SizeType i, const double val) {
            utopia_test_assert(val <= 2. + 1e-16);
        });
        
        
        double s_At_v = sum(At_v);
        utopia_test_assert(approxeq(s_At_v, size(A).get(0) * 2.));
    }
    
    
    void trilinos_apply_transpose_explicit()
    {
        auto rows = 5;
        auto cols = 6;
        TSMatrixd A = local_sparse(rows, cols, 2);
        
        {
            Write<TSMatrixd> w_A(A);
            Range r = row_range(A);
            
            for(auto i = r.begin(); i < r.end(); ++i) {
                A.set(i, i,     1.);
                A.set(i, i + 1, 1.);
            }
        }
        
        TVectord v    = local_values(rows, 1.);
        //Expilcit transpose
        TSMatrixd At  = transpose(A);
        TVectord At_v = At * v;
        
        // disp(At);
        
        each_read(At_v, [](const SizeType i, const double val) {
            utopia_test_assert(val <= 2. + 1e-16);
        });
        
        
        double s_At_v = sum(At_v);
        
        // disp(s_At_v);
        
        utopia_test_assert(approxeq(s_At_v, size(A).get(0) * 2.));
    }
    
    void trilinos_transpose()
    {
        auto rows = 5;
        auto cols = 5;
        TSMatrixd A = local_sparse(rows, cols, 2);
        auto gs = size(A);
        
        {
            Write<TSMatrixd> w_A(A);
            Range r = row_range(A);
            
            for(auto i = r.begin(); i < r.end(); ++i) {
                A.set(i, 0, 1.);
                A.set(i, gs.get(1)-1, 2.);
            }
        }
        
        auto s = size(A);
        
        TSMatrixd At = transpose(A);
        auto s_t = size(At);
        utopia_test_assert(s_t.get(0) == s.get(1));
        
        
        // disp(A);
        // std::cout << "-----------------------" << std::endl;
        // disp(At);
        // std::cout << "-----------------------" << std::endl;
        
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
        //        disp(A);
        //        disp(d);
        
        const double val = norm1(d);
        utopia_test_assert(approxeq(val, size(d).get(0)*2.-2.));
        
        TSMatrixd D = diag(d);
        TVectord x  = local_values(n, 1.);
        //        disp(D);
        //        disp(x);
        utopia_test_assert(approxeq(d, D*x));
    }
    
    void trilinos_diag_rect_matrix()
    {
        auto n = 10;
        auto m = 3;
        TSMatrixd A = local_identity(n, m);
        TVectord d;
        d  = diag(A);
        const double val = norm1(d);
        
        if( !approxeq(val, size(A).get(1) * 1.) ) {
            m_utopia_error("diag does not work on for tpetra rectangular matrices in parallel (nor serial)");
        }
    }
    
    void test_ptap(const int n, const int m)
    {
        TSMatrixd A = local_sparse(n, n, 3);
        assemble_laplacian_1D(A);
        
        TSMatrixd P;
        build_rectangular_matrix(n, m, P);
        
        TSMatrixd R_2 = transpose(P) * A;
        utopia_test_assert(R_2.implementation().is_valid(true));
        
        //For the moment this is computing (transpose(P) * A) * P
        TSMatrixd R  = utopia::ptap(A, P); //equiv: transpose(P) * A * P;
        
        utopia_test_assert(R.implementation().is_valid(true));
        
#ifdef WITH_PETSC
        //using petsc to test trilinos
        
        DSMatrixd A_petsc = local_sparse(n, n, 3);
        assemble_laplacian_1D(A_petsc);
        
        DSMatrixd P_petsc;
        build_rectangular_matrix(n, m, P_petsc);
        
        DSMatrixd R_2_petsc = transpose(P_petsc) * A_petsc;
        DSMatrixd R_petsc   = utopia::ptap(A_petsc, P_petsc);
        
        DSMatrixd R_tpetra;
        DSMatrixd R_2_tpetra;
        
        
        backend_convert_sparse(R_2, R_2_tpetra);
        backend_convert_sparse(R, R_tpetra);
        
        // disp(R_2_tpetra);
        // disp("-----------------------------");
        // disp(R_2_petsc);
        
        // write("R_t.mm", R);
        // write("R_p.m", R_petsc);
        
        // write("R2_t.mm", R_2);
        // write("R2_p.m", R_2_petsc);
        
        double diff_2 = norm2(R_2_petsc - R_2_tpetra);
        double diff   = norm2(R_petsc - R_tpetra);
        
        utopia_test_assert(approxeq(diff_2, 0.));
        utopia_test_assert(approxeq(diff, 0.));
#endif //WITH_PETSC
    }

    void test_rap(const int n, const int m)
    {
        TSMatrixd A = local_sparse(n, n, 3);
        assemble_laplacian_1D(A);
        
        TSMatrixd P;
        build_rectangular_matrix(n, m, P);
        
        TSMatrixd R = transpose(P);


        TSMatrixd res = R * A * P;

        // disp(res);
    }
    
    void trilinos_ptap_square_mat()
    {
        test_ptap(10, 10);
    }
    
    void trilinos_ptap()
    {
        //does not work in parallel
        test_ptap(10, 3);
    }

    void trilinos_rap_square_mat()
    {
        //does not work in parallel
        test_rap(10, 10);
    }

    void trilinos_rap()
    {
        //does not work in parallel
        test_rap(10, 10);
    }
    
    void trilinos_cg()
    {
        MultiLevelTestProblem<TSMatrixd, TVectord> ml_problem(10, 2);
        TVectord x = zeros(size(*ml_problem.rhs));
        (*ml_problem.rhs) *= 0.0001;
        
        ConjugateGradient<TSMatrixd, TVectord> cg;
        cg.rtol(1e-6);
        cg.atol(1e-6);
        cg.max_it(800);
        // cg.verbose(true);
        cg.update(ml_problem.matrix);

        x = *ml_problem.rhs;
        cg.apply(*ml_problem.rhs, x);
        
        double diff = norm2(*ml_problem.rhs - *ml_problem.matrix * x);;
        utopia_test_assert(approxeq(diff, 0., 1e-6));
    }
    
    template<class Matrix, class Vector>
    void test_mg()
    {
        using TransferT       = utopia::Transfer<Matrix, Vector>;
        using IPTransferT     = utopia::IPTransfer<Matrix, Vector>;
        using MatrixTransferT = utopia::MatrixTransfer<Matrix, Vector>;
        
        const static bool verbose   = false;
        const static bool use_masks = false;
        
        MultiLevelTestProblem<Matrix, Vector> ml_problem(10, 4, !use_masks);
        // ml_problem.describe();
        // ml_problem.write_matlab("./");
        
        auto smoother      = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
        auto coarse_solver = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();

        coarse_solver->set_preconditioner(std::make_shared< InvDiagPreconditioner<Matrix, Vector> >());
        coarse_solver->max_it(1000);
        
        Multigrid<Matrix, Vector> multigrid(
                                            smoother,
                                            coarse_solver
                                            );
        
        multigrid.max_it(8);
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
        
        multigrid.set_transfer_operators(transfers);
        
        Vector x = zeros(size(*ml_problem.rhs));
        multigrid.update(ml_problem.matrix);
        
        if(verbose) {
            multigrid.describe();
        }
        
        multigrid.apply(*ml_problem.rhs, x);
        
        double diff0 = norm2(*ml_problem.matrix * x);
        double diff  = norm2(*ml_problem.rhs - *ml_problem.matrix * x);
        double rel_diff = diff/diff0;
        
        utopia_test_assert(rel_diff < 1e-8);
    }
    
    void trilinos_local_row_view()
    {
        auto rows = 3;
        auto cols = 4;
        TSMatrixd A = local_sparse(rows, cols, 2);
        
        {
            Write<TSMatrixd> w_A(A);
            Range r = row_range(A);
            
            for(auto i = r.begin(); i < r.end(); ++i) {
                A.set(i, i,     i);
                A.set(i, i + 1, i + 1);
            }
        }
        
        TSMatrixd At = transpose(A);
        
        auto &M = A;
        
        auto rr = row_range(M);
        for(auto i = rr.begin(); i < rr.end(); ++i) {
            RowView<TSMatrixd> row(M, i, true);
            for(auto j = 0; j < row.n_values(); ++j) {
                int col = row.col(j);
                int val = row.get(j);
                utopia_test_assert(col == val);
            }
        }
    }
    
    void trilinos_range()
    {
        SizeType n = 10, m = 5;
        SizeType gn = mpi_world_size() * n, gm = mpi_world_size() * m;
        
        TSMatrixd P;
        build_rectangular_matrix(n, m, P);
        
        Range rr = row_range(P);
        
        utopia_test_assert(rr.begin() == n * mpi_world_rank());
        utopia_test_assert(rr.end()   == n * (mpi_world_rank() + 1));
        utopia_test_assert(gn == size(P).get(0));
        utopia_test_assert(gm == size(P).get(1));
    }
    
    void trilinos_e_mul()
    {
        
        int n = 10;
        TVectord v    = local_values(n, 1.);
        TVectord ones = local_values(local_size(v).get(0), 1.);
        
        
        TVectord ones_mul_v = e_mul(ones, v);
        v = e_mul(ones, v);
        
        double sv = sum(v);
        utopia_test_assert(approxeq(sv, n * mpi_world_size()));
    }


    template<class Matrix, class Vector>
    void st_cg_test()
    {
        typename Vector::SizeType _n = 10; 

        SteihaugToint<Matrix, Vector, HOMEMADE> cg;
        cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());
        // cg.set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());

        cg.rtol(1e-7);
        cg.atol(1e-6);
        cg.max_it(_n);
        cg.verbose(false);

        Matrix A = sparse(_n, _n, 3);
        assemble_symmetric_laplacian_1D(A, true);

        Vector rhs = values(_n, 975.9);

        {
            auto r = range(rhs);

            Write<Vector> w(rhs);
            if(r.begin() == 0)  rhs.set(0, 0.0); 
            if(r.end()   == _n) rhs.set(_n-1, 0.0); 
        }           

        Vector x = zeros(size(rhs));

        cg.tr_constrained_solve(A, -1.0 * rhs, x, 1e15);
        utopia_test_assert(approxeq(rhs, A * x, 1e-5));

    }

    void stcg_pt_test()
    {

#ifdef WITH_PETSC
        //petsc version
        st_cg_test<DSMatrixd, DVectord>();
#endif //WITH_PETSC
        st_cg_test<TSMatrixd, TVectord>();
    }


    void trilinos_mg_1D()
    {
        // if(mpi_world_size() > 1) return;
        //petsc version
#ifdef WITH_PETSC
        test_mg<DSMatrixd, DVectord>();
#endif //WITH_PETSC
        //trilinos version
        test_mg<TSMatrixd, TVectord>();
    }
    
    
    void trilinos_mg()
    {
        // if(mpi_world_size() > 1) return;
        
        using MatrixT = utopia::TSMatrixd;
        using VectorT = utopia::TVectord;
        
        // using MatrixT = utopia::DSMatrixd;
        // using VectorT = utopia::DVectord;
        
        VectorT rhs;
        MatrixT A, I;
        
        Multigrid<MatrixT, VectorT> multigrid(
                                              std::make_shared<ConjugateGradient<MatrixT, VectorT, HOMEMADE>>(),
                                              std::make_shared<ConjugateGradient<MatrixT, VectorT, HOMEMADE>>()
                                              // std::make_shared<SOR<MatrixT, VectorT>>(),
                                              // std::make_shared<Factorization<MatrixT, VectorT>>()
                                              );
        
#ifdef WITH_PETSC
        
        bool ok = true;
        //FIXME needs trilinos formats but for the moment lets use petsc's
        {
            DSMatrixd petsc_A, petsc_I;
            DVectord petsc_rhs;
            
            const std::string folder =  Utopia::instance().get("data_path") + "/laplace/matrices_for_petsc";
            
            ok = read(folder + "/f_rhs", petsc_rhs); utopia_test_assert(ok);
            ok = read(folder + "/f_A", petsc_A);     utopia_test_assert(ok);
            ok = read(folder + "/I_3", petsc_I);     utopia_test_assert(ok);
            
            backend_convert_sparse(petsc_I, I);
            backend_convert_sparse(petsc_A, A);
            backend_convert(petsc_rhs, rhs);
        }
        
        // write("A.mm", A);
        // write("I.mm", I);
        
        std::vector<std::shared_ptr<MatrixT>> interpolation_operators;
        interpolation_operators.push_back(make_ref(I));
        
        multigrid.set_transfer_operators(std::move(interpolation_operators));
        multigrid.max_it(20);
        multigrid.atol(1e-15);
        multigrid.stol(1e-15);
        multigrid.rtol(1e-15);
        // multigrid.verbose(true);
        multigrid.set_fix_semidefinite_operators(true);
        multigrid.must_generate_masks(true);
        VectorT x = local_zeros(local_size(rhs));
        
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
            // auto val = row.get(0);
            
            utopia_test_assert(col == i || col == i - 1 || col == i  + 1);
        }
        
        
        TSMatrixd B = local_sparse(4, 4, 3);
        
        {
            Write<TSMatrixd> w_(B, GLOBAL_ADD);
            
            auto r = row_range(B);
            auto s = size(B);
            
            for(auto i = r.begin(); i < r.end(); ++i) {
                B.set(i, i, i);
                
                if(i + 1 < s.get(1)) {
                    B.set(i, i+1, i+1);
                }
                
                if(i - 1 >= 0) {
                    B.set(i, i-1, i-1);
                }
            }
        }
        
        each_read(B, [](const SizeType i, const SizeType j, const double val) {
            SizeType j_val = val;
            utopia_test_assert(j_val == val);
        });
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
    
    void trilinos_exp()
    {
#ifdef WITH_PETSC
        TVectord x = local_values(10, 2.);
        DVectord y = local_values(10, 2.);
        
        TVectord ex = exp(x);
        DVectord ey = exp(y);
        
        utopia_test_assert(cross_backend_approxeq(ey, ex));
        
#endif //WITH_PETSC
    }
    
    void trilinos_diag_ops()
    {
        int n = 10;
        
        TSMatrixd m = local_sparse(n, n, 3);
        assemble_laplacian_1D(m);
        
        TSMatrixd m_copy = m;
        
        TVectord d  = local_values(n, -1.);
        TSMatrixd D = diag(d);
        m += D;
        
        TVectord ones = local_values(n, 1);
        
        double sum_m      = sum(m * ones);
        double sum_d      = sum(d);
        double sum_D      = sum(D * ones);
        double sum_m_copy = sum(m_copy * ones);
        
        utopia_test_assert(approxeq(sum_d, sum_D));
        utopia_test_assert(approxeq(sum_m, sum_d + sum_m_copy));
        
        TSMatrixd m_new = m_copy + D;
        double sum_m_new = sum(m * ones);
        
        utopia_test_assert(approxeq(sum_m_new, sum_d + sum_m_copy));
    }
    
    void trilinos_bratu_1D()
    {
#ifdef WITH_PETSC
        int n = 10;
        double lambda = 1.9;
        
        auto fun_tpetra = std::make_shared<Bratu1D<TSMatrixd, TVectord> >(n, lambda);
        auto fun_petsc  = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(n, lambda);
        
        auto init_expr = values(n, 1.);
        
        TVectord x_tpetra = init_expr;
        DVectord x_petsc  = init_expr;
        
        double val_tpetra = 0.;
        double val_petsc  = 0.;
        
        bool ok = true;
        ok = fun_tpetra->value(x_tpetra, val_tpetra); assert(ok);
        ok = fun_petsc->value(x_petsc,   val_petsc);  assert(ok);
        
        utopia_test_assert(cross_backend_approxeq(x_petsc, x_tpetra));
        utopia_test_assert(approxeq(val_tpetra, val_petsc));
        
        TVectord grad_tpetra;
        DVectord grad_petsc;
        
        ok = fun_tpetra->gradient_no_rhs(x_tpetra, grad_tpetra); assert(ok);
        ok = fun_petsc->gradient_no_rhs(x_petsc, grad_petsc);    assert(ok);
        
        utopia_test_assert(cross_backend_approxeq(grad_petsc, grad_tpetra));
        
        //last part fails
        
        TSMatrixd H_tpetra;
        DSMatrixd H_petsc;
        
        ok = fun_tpetra->hessian(x_tpetra, H_tpetra); assert(ok);
        ok = fun_petsc->hessian(x_petsc, H_petsc);    assert(ok);
        
        // write("H_p.m", H_petsc);
        // write("H_t.m", H_tpetra);
        
        
        DSMatrixd H_converted;
        backend_convert_sparse(H_tpetra, H_converted);
        
        // write("H_c.m", H_converted);
        
        utopia_test_assert(cross_backend_approxeq(H_petsc, H_tpetra));
        
#endif //WITH_PETSC
    }
    
    template<class Matrix, class Vector>
    void rmtr_test()
    {
        using IPTransferT = utopia::IPTransfer<Matrix, Vector>;
        
        BratuMultilevelTestProblem<Matrix, Vector> problem(2, true);
        problem.verbose = false;
        // problem.verbose = tre;
        
        Vector x = values(problem.n_dofs[problem.n_levels -1 ], 0.0);
        
        std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >  level_functions(problem.n_levels);
        
        for(auto l = 0; l < problem.n_levels; l++)
        {
            auto fun = std::make_shared<Bratu1D<Matrix, Vector> >(problem.n_dofs[l]);
            level_functions[l] = fun;
            
            // making sure that fine level IG is feasible
            if(l + 1 == problem.n_levels) {
                fun->apply_bc_to_initial_guess(x);
            }
        }
        
        auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
        tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());

        tr_strategy_coarse->atol(1e-12);
        tr_strategy_coarse->rtol(1e-12);
        
        auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();

        
        tr_strategy_fine->set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());
        tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());

        tr_strategy_fine->atol(1e-12);
        tr_strategy_fine->rtol(1e-12);
        
        // auto rmtr = std::make_shared<RMTR<Matrix, Vector, SECOND_ORDER>  >(tr_strategy_coarse, tr_strategy_fine);
        auto rmtr = std::make_shared<RMTR<Matrix, Vector, GALERKIN> >(tr_strategy_coarse, tr_strategy_fine);
        std::vector< std::shared_ptr<Transfer<Matrix, Vector>> > transfers;
        for(std::size_t i = 0; i < problem.prolongations.size(); ++i) {
            transfers.push_back( std::make_shared<IPTransferT>( problem.prolongations[i], 0.5) );
        }
        
        rmtr->set_transfer_operators(transfers);
        
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
        
        
        rmtr->handle_equality_constraints();
        bool ok = rmtr->solve(x);
        
        utopia_test_assert(ok);
    }

    void trilinos_replace_value()
    {
        TSMatrixd A = local_sparse(2, 3, 4);

        auto gs = size(A);
        auto r = row_range(A);

        {
            Write<TSMatrixd> w_(A);

            for(auto i = r.begin(); i < r.end(); ++i) {
                // A.set(i, i, 1.);
                A.set(i, i + 1, 2.);
                A.set(i, 0, 3.);
                A.set(i, std::max(gs.get(1)-i-1, SizeType(0)), 4.);
            }
        }

        {
            Write<TSMatrixd> w_(A);

            for(auto i = r.begin(); i < r.end(); ++i) {
                A.set(i, 0, -3.);
                A.set(i, std::max(gs.get(1)-i-1, SizeType(0)), -4.);
            }
        }
    }
    
    
    void trilinos_copy_write()
    {   
        TSMatrixd P;
        build_rectangular_matrix(10, 5, P);

        P = transpose(P);
        
        TSMatrixd P2 = P;
        P2 *= 0.;

        {
            Write<TSMatrixd> w_(P2);
            each_read(P, [&P2](const SizeType i, const SizeType j, const double value) {
                P2.set(i, j, value * 2.);
            });
        }

        TSMatrixd A;
        MultiLevelTestProblem<TSMatrixd, TVectord> ml_problem(10, 2);
        A = *ml_problem.interpolators[0];

        A = transpose(A);

        TSMatrixd A2 = A;
        
        A2 *= 0.;

        {
            Write<TSMatrixd> w_(A2);
            each_read(A, [&A2](const SizeType i, const SizeType j, const double value) {
                A2.set(i, j, 1.);
            });
        }
    }

#ifdef WITH_PETSC
    void trilinos_copy_write_big()
    {   
        DSMatrixd petsc_P;

        const std::string folder =  Utopia::instance().get("data_path") + "/laplace/matrices_for_petsc";
        bool ok = read(folder + "/I_2", petsc_P); utopia_test_assert(ok);

        TSMatrixd P;
        backend_convert_sparse(petsc_P, P);

        P = transpose(P);
        
        TSMatrixd P2 = P;
        P2 *= 0.;

        {
            Write<TSMatrixd> w_(P2);
            each_read(P, [&P2](const SizeType i, const SizeType j, const double value) {
                P2.set(i, j, value * 2.);
            });
        }
    }
#endif //WITH_PETSC
    
    void trilinos_ghosted()
    {
        const int n   = mpi_world_size() * 2;
        const int off = mpi_world_rank() * 2;
        
        std::vector<TVectord::SizeType> ghosts{ (off + 3) % n };
        TVectord v = ghosted(2, n, ghosts);
        
        auto r = range(v);
        
        {
            Write<TVectord> w_v(v);
            for(auto i = r.begin(); i != r.end(); ++i) {
                v.set(i, i);
            }
        }
        
        // synchronize(v);
        
        {
            Read<TVectord> r_v(v);
            std::vector<SizeType> index{(off + 3) % n};
            std::vector<double> values;
            v.get(index, values);
            utopia_test_assert(index[0] == SizeType(values[0]));
        }
        
        // disp(v);
    }
    
    void trilinos_rmtr()
    {
#ifdef WITH_PETSC
        //petsc version
        rmtr_test<DSMatrixd, DVectord>();
#endif //WITH_PETSC
        
        rmtr_test<TSMatrixd, TVectord>();
    }
    
    void trilinos_matrix_norm()
    {
        TSMatrixd m = local_identity(10, 10);
        
        double nm = norm2(m);
        utopia_test_assert( approxeq(nm, std::sqrt(1.*size(m).get(0))) );
    }

#ifdef HAVE_BELOS_TPETRA

    void trilinos_belos()
    {
        std::string xml_file = Utopia::instance().get("data_path") + "/UTOPIA_belos.xml";
        
        Parameters params;
        params.set_param_file_name(xml_file);
        BelosSolver<TSMatrixd, TVectord> solver(params);

        MultiLevelTestProblem<TSMatrixd, TVectord> ml_problem(10, 2);
        TVectord x = zeros(size(*ml_problem.rhs));
        (*ml_problem.rhs) *= 0.0001;
        
        double diff0 = norm2(*ml_problem.rhs - *ml_problem.matrix * x);

        x = *ml_problem.rhs;
        solver.solve(*ml_problem.matrix, *ml_problem.rhs, x);
        
        double diff = norm2(*ml_problem.rhs - *ml_problem.matrix * x);

        utopia_test_assert(approxeq(diff/diff0, 0., 1e-6));
    }

#endif //HAVE_BELOS_TPETRA

#ifdef HAVE_AMESOS2_KOKKOS

    void trilinos_amesos2()
    {
        std::string xml_file = Utopia::instance().get("data_path") + "/UTOPIA_amesos.xml";
        
        Parameters params;
        params.set_param_file_name(xml_file);
        Amesos2Solver<TSMatrixd, TVectord> solver(params);

        MultiLevelTestProblem<TSMatrixd, TVectord> ml_problem(10, 2);
        TVectord x = zeros(size(*ml_problem.rhs));
        (*ml_problem.rhs) *= 0.0001;
        
        double diff0 = norm2(*ml_problem.rhs - *ml_problem.matrix * x);

        solver.solve(*ml_problem.matrix, *ml_problem.rhs, x);
        
        double diff  = norm2(*ml_problem.rhs - *ml_problem.matrix * x);

        utopia_test_assert(approxeq(diff/diff0, 0., 1e-6));
    }

#endif //HAVE_AMESOS2_KOKKOS



#ifdef WITH_PETSC
    void trilinos_transform()
    {
        DSMatrixd petsc_P;

        const std::string folder =  Utopia::instance().get("data_path") + "/laplace/matrices_for_petsc";
        bool ok = read(folder + "/I_2", petsc_P); utopia_test_assert(ok);

        TSMatrixd P;
        backend_convert_sparse(petsc_P, P);

        double sum_P = sum(P);
        P = transpose(P);

        each_apply(P, [](const double value) -> double {
            return value * 2.;
        });

        double sum_P_2 = sum(P);

        utopia_test_assert(approxeq(sum_P * 2., sum_P_2));


        each_transform(P, [](const SizeType i, const SizeType j, const double) -> double {
            return j;
        });

        each_read(P, [](const SizeType i, const SizeType j, const double value) {
            if(j != SizeType(value)) {
                std::cout << i << " " << j << " " << value << std::endl;
            }

            utopia_test_assert(j == SizeType(value));
        });
    }
#endif //WITH_PETSC

    void trilinos_set_zeros()
    {
        TSMatrixd m = local_identity(10, 10);
        using SizeT = UTOPIA_SIZE_TYPE(TSMatrixd);
        auto rr = row_range(m);

        std::vector<SizeT> index;
        index.push_back(rr.begin());
        set_zero_rows(m, index, 2.);

        // disp(m);

        Read<TSMatrixd> r(m);
        auto val = m.get(rr.begin(), rr.begin());
        utopia_test_assert(val == 2.);
    }
    
    void run_trilinos_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("TrilinosTest");

        UTOPIA_RUN_TEST(stcg_pt_test); 

        UTOPIA_RUN_TEST(trilinos_structure);
        UTOPIA_RUN_TEST(trilinos_build);
        UTOPIA_RUN_TEST(trilinos_build_identity);
        UTOPIA_RUN_TEST(trilinos_accessors);
        UTOPIA_RUN_TEST(trilinos_vec_scale);
        UTOPIA_RUN_TEST(trilinos_mat_scale);
        UTOPIA_RUN_TEST(trilinos_vec_axpy);
        UTOPIA_RUN_TEST(trilinos_vec_minus);
        UTOPIA_RUN_TEST(trilinos_mv);
        UTOPIA_RUN_TEST(trilinos_mm);
        UTOPIA_RUN_TEST(trilinos_m_tm);
        UTOPIA_RUN_TEST(trilinos_diag);
        UTOPIA_RUN_TEST(trilinos_read);
        
        UTOPIA_RUN_TEST(trilinos_rect_matrix);
        UTOPIA_RUN_TEST(trilinos_e_mul);
        UTOPIA_RUN_TEST(trilinos_row_view);
        UTOPIA_RUN_TEST(trilinos_apply_transpose);
        UTOPIA_RUN_TEST(trilinos_set);
        UTOPIA_RUN_TEST(trilinos_residual);
        
        UTOPIA_RUN_TEST(trilinos_matrix_access);
        UTOPIA_RUN_TEST(trilinos_matrix_norm);
        UTOPIA_RUN_TEST(trilinos_diag_ops);
        
        UTOPIA_RUN_TEST(trilinos_rmtr);
        
        UTOPIA_RUN_TEST(trilinos_transpose);
        
        UTOPIA_RUN_TEST(trilinos_apply_transpose_explicit);
        UTOPIA_RUN_TEST(trilinos_each_read_transpose);
        UTOPIA_RUN_TEST(trilinos_local_row_view);
        UTOPIA_RUN_TEST(trilinos_ptap_square_mat);
        
        UTOPIA_RUN_TEST(trilinos_range);
        UTOPIA_RUN_TEST(trilinos_mg_1D);
        UTOPIA_RUN_TEST(trilinos_bratu_1D);
        UTOPIA_RUN_TEST(trilinos_replace_value);
        UTOPIA_RUN_TEST(trilinos_ghosted);  
        UTOPIA_RUN_TEST(trilinos_set_zeros);
        UTOPIA_RUN_TEST(trilinos_copy_write);

        ////////////////////////////////////////////
        //test that fail on GPU if the env variables are not set correctly for cuda
        UTOPIA_RUN_TEST(trilinos_exp);
        UTOPIA_RUN_TEST(trilinos_mg);
        UTOPIA_RUN_TEST(trilinos_cg);
        UTOPIA_RUN_TEST(trilinos_ptap);
        ////////////////////////////////////////////

        UTOPIA_RUN_TEST(trilinos_rap);
        UTOPIA_RUN_TEST(trilinos_rap_square_mat);
        

#ifdef HAVE_BELOS_TPETRA
        UTOPIA_RUN_TEST(trilinos_belos);
#endif //HAVE_BELOS_TPETRA  

//#ifdef HAVE_AMESOS2_TPETRA
        UTOPIA_RUN_TEST(trilinos_amesos2);
//#endif //HAVE_AMESOS2_TPETRA

#ifdef WITH_PETSC
        UTOPIA_RUN_TEST(trilinos_transform);
        UTOPIA_RUN_TEST(trilinos_petsc_interop);
        UTOPIA_RUN_TEST(trilinos_copy_write_big);
#endif //WITH_PETSC    

        //Fails on multinode GPU
        UTOPIA_RUN_TEST(trilinos_mat_axpy);
        ////////////////////////////////////////////

        if(mpi_world_size() <= 3) {
            //working up to 3 processes
            UTOPIA_RUN_TEST(trilinos_row_view_and_loops);
        }

        //tests that fail in parallel
        // if(mpi_world_size() == 1) {
               
        // } else {
        //     m_utopia_warning_once("several tests left out for parallel execution");
        // }
        
        //tests that always fail
        // UTOPIA_RUN_TEST(trilinos_diag_rect_matrix);
        
        UTOPIA_UNIT_TEST_END("TrilinosTest");
    }
}

#else //WITH_TRILINOS

namespace utopia
{
    void run_trilinos_test() {}
}

#endif //WITH_TRILINOS


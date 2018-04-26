#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS

#include "utopia_TrilinosTest.hpp"
#include "utopia.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_solvers.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
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
        assert(approxeq(nz, 0.));
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
        
        assert(s.get(0) == n * mpi_world_size());
        assert(ls.get(0) == n);

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
        assert(approxeq(val, n * mpi_world_size() * 1.5));
    }

    void trilinos_vec_scale()
    {
        auto n = 10;
        TVectord y = local_values(n, 1.);

        y *= 2.;

        double val = norm1(y);
        assert(approxeq(val, size(y).get(0) * 2.));


        TVectord y2 = y * 2.;

        val = norm1(y2);
        assert(approxeq(val, size(y2).get(0) * 4.));
    }

    void trilinos_mat_scale()
    {
        auto n = 10;
        TSMatrixd Y = local_sparse(n, n, 3);
        assemble_laplacian_1D(Y);

        Y *= 2.;

        TVectord x = local_values(n, 1.);
        double val = norm1(Y * x);
        assert(approxeq(val, 0.));


        TSMatrixd Y2 = Y * 2.;

        val = norm1(Y2 * x);
        assert(approxeq(val, 0.));
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
        assert(approxeq(val, 0.));
    }

    void trilinos_mv()
    {
        auto n = 10;
        TVectord x = local_values(n, 5.);
        TSMatrixd m = sparse(n, n, 3);
        assemble_laplacian_1D(m);

        TVectord y = m * x;

        const double val = norm2(y);
        assert(approxeq(val, 0.));
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
        // disp(A);
        // disp(At);

        auto s_t = size(At);
        assert(s_t.get(0) == s.get(1));
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
        assert(approxeq(val, 0.));
    }

    template<class Matrix>
    static void build_rectangular_matrix(const SizeType &n, const SizeType &m, Matrix &mat)
    {
        mat  = local_sparse(n, m, 1);

        Write<TSMatrixd> w_(mat);
        auto r = row_range(mat);
        auto cols = size(mat).get(1);
        for(auto i = r.begin(); i < r.end(); ++i) {
            // if(i >= cols) {
            //     break;
            // }

            if(i < cols) {
                mat.set(i, i, 1.);
            } else {
                mat.set(i, 0, 1e-16);
            }
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
        assert(approxeq(val, size(d).get(0)*2.-2.));

        TSMatrixd D = diag(d);
        TVectord x  = local_values(n, 1.);
        assert(approxeq(d, D*x));
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
        TSMatrixd R = utopia::ptap(A, P);

        // disp(A);
        // disp(P);
        // disp(R);

        //FIXME write test here
    }

    void trilinos_mg()
    {
        if(mpi_world_size() > 1) return;

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
        
            const std::string folder =  Utopia::instance().get("data_path") + "/laplace/matrices_for_petsc";
            
            ok = read(folder + "/f_rhs", petsc_rhs); assert(ok);
            ok = read(folder + "/f_A", petsc_A);     assert(ok);
            // ok = read(folder + "/I_2", I_2);   assert(ok);
            ok = read(folder + "/I_3", petsc_I);     assert(ok);
           
            backend_convert_sparse(petsc_I, I);
            backend_convert_sparse(petsc_A, A);
            backend_convert(petsc_rhs, rhs);
        }

        std::vector<std::shared_ptr<TSMatrixd>> interpolation_operators;
        interpolation_operators.push_back(make_ref(I));
        
        multigrid.set_transfer_operators(std::move(interpolation_operators));
        multigrid.max_it(20);
        multigrid.atol(1e-15);
        multigrid.stol(1e-15);
        multigrid.rtol(1e-15);
        // multigrid.verbose(true);
        multigrid.must_generate_masks(false);
        TVectord x = local_zeros(local_size(rhs));

        try {
            multigrid.update(make_ref(A));
            ok = multigrid.apply(rhs, x); assert(ok);
        } catch(const std::exception &ex) {
            std::cout << ex.what() << std::endl;
            assert(false);
        }

        std::cout << std::flush;
        assert(approxeq(rhs, A * x, 1e-6));
#endif //WITH_PETSC
        
    }

    void row_view_and_loops()
    {
        auto n = 10;
        auto m = 3;

        TSMatrixd P;
        build_rectangular_matrix(n, m, P);

        each_read(P, [](const SizeType i, const SizeType j, const double val) {

        });

        TSMatrixd P_t = transpose(P);


        each_read(P_t, [](const SizeType i, const SizeType j, const double val) {

        });
    }

    void trilinos_read()
    {
        TSMatrixd m;
        auto path = Utopia::instance().get("data_path") + "/matrixmarket/gre_343_343_crg.mm";
        bool ok = read(path, m);
        assert(ok);
    }
    
    void run_trilinos_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("TrilinosTest");
        UTOPIA_RUN_TEST(trilinos_build);
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
        UTOPIA_RUN_TEST(trilinos_mg);

        //does not work
        // UTOPIA_RUN_TEST(row_view_and_loops);
        UTOPIA_UNIT_TEST_END("TrilinosTest");
    }
}

#else //WITH_TRILINOS

namespace utopia
{
    void run_trilinos_test() {}
}

#endif //WITH_TRILINOS

#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS

#include "utopia_TrilinosTest.hpp"
#include "utopia.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_solvers.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

namespace utopia {   
 
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

    void trilinos_mg()
    {
        if(mpi_world_size() > 1) return;

        TVectord rhs;
        TSMatrixd A, I;

        Multigrid<TSMatrixd, TVectord> multigrid(
            std::make_shared<ConjugateGradient<TSMatrixd, TVectord>>(),
            std::make_shared<ConjugateGradient<TSMatrixd, TVectord>>()
        );
        
        const std::string data_path = Utopia::instance().get("data_path");
        const std::string folder = data_path + "/mg";
        
        //needs trilinos formats
        // read(folder + "/rhs.bin", rhs);
        // read(folder + "/A.bin", A);
        // read(folder + "/I.bin", I);
        
        std::vector<std::shared_ptr<TSMatrixd>> interpolation_operators;
        interpolation_operators.push_back(make_ref(I));
        
        multigrid.init_transfer_from_fine_to_coarse(std::move(interpolation_operators));
        multigrid.max_it(20);
        multigrid.atol(1e-15);
        multigrid.stol(1e-15);
        multigrid.rtol(1e-15);
        // multigrid.verbose(verbose);
        
        TVectord x = zeros(A.size().get(0));

        int block_size = 2;
        multigrid.block_size(block_size);
        multigrid.update(make_ref(A));
        multigrid.apply(rhs, x);

        assert(approxeq(rhs, A * x, 1e-6));
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
        UTOPIA_RUN_TEST(trilinos_vec_axpy);
        UTOPIA_RUN_TEST(trilinos_mat_axpy);
        UTOPIA_RUN_TEST(trilinos_mv);
        UTOPIA_RUN_TEST(trilinos_mm);
        UTOPIA_RUN_TEST(trilinos_diag);
        UTOPIA_RUN_TEST(trilinos_read);
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

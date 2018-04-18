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
        
        // double val = norm1(y);
        // assert(approxeq(val, n * mpi_world_size() * 1.5));
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
    
    void run_trilinos_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("TrilinosTest");
        UTOPIA_RUN_TEST(trilinos_build);
        UTOPIA_RUN_TEST(trilinos_accessors);
        UTOPIA_RUN_TEST(trilinos_vec_axpy);
        UTOPIA_RUN_TEST(trilinos_mat_axpy);
        UTOPIA_RUN_TEST(trilinos_mv);
        UTOPIA_UNIT_TEST_END("TrilinosTest");
    }
}

#else //WITH_TRILINOS

namespace utopia
{
    void run_trilinos_test() {}
}

#endif //WITH_TRILINOS

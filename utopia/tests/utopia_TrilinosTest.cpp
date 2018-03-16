#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS

#include "utopia_TrilinosTest.hpp"
#include "utopia.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_solvers.hpp"

namespace utopia {   

    template<class Matrix>
    static void build_laplacian(const int n, Matrix &m)
    {
        m = local_sparse(n, n, 3);

        Write<TSMatrixd> w_m(m);
        Range r = row_range(m);
        Size s = size(m);

        for(auto i = r.begin(); i < r.end(); ++i) {
            if(i > 0) {
                m.set(i, i - 1, -1.);
            }  

            if(i + 1 < s.get(1)) {
                m.set(i, i + 1, -1.);
            }

            if(i == 0 || i == n - 1) {
                m.set(i, i, 1);
            } else {
                m.set(i, i, 2.);
            }
        }
    }
    
    void trilinos_build_test()
    {
        auto n = 10;
        TVectord v = local_values(n, 1.);
        
        //FIXME replace this with an actual test
        // disp(v);

        TSMatrixd m;
        build_laplacian(n, m);

        //FIXME replace this with an actual test
        // disp(m);
    }
    
    void trilinos_accessors_test()
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
    
    void trilinos_axpy_test()
    {
        auto n = 10;
        TVectord y = local_values(n, 1.);
        TVectord x = local_values(n, 5.);
        auto alpha = 0.1;
        y += alpha * x;
        
        double val = norm1(y);
        assert(approxeq(val, n * mpi_world_size() * 1.5));
        
        //FIXME replace this with an actual test
        // disp(y);
    }

    void trilinos_mult_test()
    {
        auto n = 10;
        TVectord x = local_values(n, 5.);
        TSMatrixd m;
        build_laplacian(n, m);

        TVectord y = m * x;

        const double val = norm2(y);
        assert(approxeq(val, 0.));
    }
    
    void run_trilinos_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("TrilinosTest");
        UTOPIA_RUN_TEST(trilinos_build_test);
        UTOPIA_RUN_TEST(trilinos_accessors_test);
        UTOPIA_RUN_TEST(trilinos_axpy_test);
        UTOPIA_RUN_TEST(trilinos_mult_test);
        UTOPIA_UNIT_TEST_END("TrilinosTest");
    }
}

#else //WITH_TRILINOS

namespace utopia
{
    void run_trilinos_test() {}
}

#endif //WITH_TRILINOS

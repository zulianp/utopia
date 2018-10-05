
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS

#include "utopia_KokkosTest.hpp"
#include "utopia.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_solvers.hpp"
#include "utopia_Each_Parallel.hpp"

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

namespace utopia {


    void kokkos_max()
    {
        auto n = 10;
        TVectord v = local_values(n, 1.);
        TVectord w = local_values(n, 10.);

        TVectord z = max(w , v);
        
        disp("max_between vectors");
        disp(z);

    }
    
    void kokkos_min()
    {
        auto n = 10;
        TVectord v = local_values(n, 1.);
        TVectord w = local_values(n, 10.);
        
        TVectord z = min(w , v);
        
        disp("min_between vectors");
        disp(z);
        
    }


    void kokkos_sum_reduction()
    {
        auto n = 10;
        
        TVectord v = local_values(n, 1.);
        
        double z = sum(v);
        
        disp("sum");
        disp(z);
        
    }
    
    void kokkos_min_reduction()
    {
        auto n = 10;
        TVectord v = local_values(n, 1.);
        auto r = range(v);
       
        each_write(v, [&r](const SizeType i) -> double {
            return i - r.begin();
        });
        
        double z = min(v);
        disp("min");
        disp(z);
        
    }
    
    void kokkos_max_reduction()
    {
        auto n = 10;
        TVectord v = local_values(n, 1.);
        auto r = range(v);
        
        each_write(v, [&r](const SizeType i) -> double {
            return i - r.begin();
        });
//        disp(v);
        double z = max(v);
        disp("max");
        disp(z);
        
    }

       void kokkos_write(){
    
        // const int N = 10;
        // typedef Kokkos::View<double*[1]> view_type;
    
        // view_type a ("A", N);
        // int b = 10;
        // KokkosWrite<view_type ,double>w(a,b);
        
        // Kokkos::parallel_for (N, w);
        //std::cout<<"I am writing TEST"<<'\n';
        auto n=10;

        TVectord w = local_values(n, 10);

        auto r = range(w);

        each_write_parallel(w, KOKKOS_LAMBDA(const SizeType i) -> double {
            return i ;
        });

        std::cout<<"Disp vector"<<'\n';

        disp(w);
    }


        void kokkos_read(){
    
        // const int N = 10;
        // typedef Kokkos::View<double*[1]> view_type;
    
        // view_type a ("A", N);
        // int b = 10;
        // KokkosWrite<view_type ,double>w(a,b);
        
        // Kokkos::parallel_for (N, w);
        //std::cout<<"I am reading TEST"<<'\n';

        auto n=10;

        TVectord w = local_values(n, 50);
        auto r = range(w);
        each_read_parallel(w, KOKKOS_LAMBDA(const SizeType i, const double entry) 
            { });
    }

  
    void run_kokkos_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("KokkosTest");
        UTOPIA_RUN_TEST(kokkos_max);
        UTOPIA_RUN_TEST(kokkos_min);
        UTOPIA_RUN_TEST(kokkos_write);
        UTOPIA_RUN_TEST(kokkos_read);
        UTOPIA_RUN_TEST(kokkos_min_reduction);
        UTOPIA_RUN_TEST(kokkos_max_reduction);
        UTOPIA_RUN_TEST(kokkos_sum_reduction);
        UTOPIA_UNIT_TEST_END("KokkosTest");
    }
}

#else //WITH_TRILINOS

namespace utopia
{
    void run_kokkos_test() {}
}

#endif //WITH_TRILINOS

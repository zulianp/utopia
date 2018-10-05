
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

#include "utopia_Trilinos_Each_Parallel.hpp"
#include <cmath>

namespace utopia {


    void kokkos_max()
    {
        auto n = 10;

        TVectord v = local_values(n, 1.);

        TVectord w = local_values(n, 10.);

        TVectord z = max(w , v);
        
        disp("kokkos_min_binary");

        disp(z);

    }
    
    void kokkos_min()
    {
        auto n = 10;

        TVectord v = local_values(n, 1.);

        TVectord w = local_values(n, 10.);
        
        TVectord z = min(w , v);
        
        disp("kokkos_min_binary");

        disp(z);
        
    }


    void kokkos_sum_reduction()
    {
        auto n = 10;
        
        TVectord v = local_values(n, 1.);
        
        double z = sum(v);
        
        disp("kokkos_sum_reduction");
        
        disp(z);
        
    }
    
    void kokkos_min_reduction()
    {
        auto n = 10;

        TVectord v = local_values(n, 1.);

        auto r = range(v);
       
        each_write(v,  KOKKOS_LAMBDA(const SizeType i) -> double {
            return i - r.begin();
        });
        
        double z = min(v);

        disp("kokkos_min_reduction");

        disp(z);
        
    }
    
    void kokkos_max_reduction()
    {
        auto n = 10;

        TVectord v = local_values(n, 1.);

        auto r = range(v);
        
        each_write(v, KOKKOS_LAMBDA(const SizeType i) -> double {
            return i - r.begin();
        });

        double z = max(v);

        disp("kokkos_max_reduction");

        disp(z);
        
    }


    void kokkos_write(){
    
        auto n=10;

        TVectord w = local_values(n, 10);

        each_write_parallel(w, KOKKOS_LAMBDA(const SizeType i) -> double {
            return i ;
        });

        disp("kokkos_write");

        disp(w);
    }


    void kokkos_read(){

        auto n=10;

        TVectord w = local_values(n, 50);

        each_read_parallel(w, KOKKOS_LAMBDA(const SizeType i, const double entry) 
            { });
    }


    void kokkos_apply()
    {
        auto nr = 3;
        auto nc = 3;

        TSMatrixd P = local_sparse(nr, nc, 1);

    {
        Write<TSMatrixd> w_(P);
        auto r = row_range(P);
        auto cols = size(P).get(1);
        for(auto i = r.begin(); i < r.end(); ++i) {
            if(i >= cols) {
                break;
            }
            
            P.set(i, i, 1.);
            
        }
    }

        disp(P.size().get(0));


        auto P_trilinos = raw_type(P);

        auto n_rows_trilinos=P_trilinos->getGlobalNumRows();

        std::cout<<"size_rows_trilinos==> "<<n_rows_trilinos<<std::endl;

        auto local_mat = P_trilinos->getLocalMatrix();

        auto n = local_mat.numRows();

        disp("size_rows_Kokkos");
        
        disp(n);

        each_apply_parallel(P, [](const double value) -> double {
            return value * 2.;
        });


        Kokkos::parallel_for( team_policy( n, Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type &teamMember) {
            
            const int j = teamMember.league_rank();
            auto row = local_mat.row(j);
            auto n_values = row.length;

            Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember, n_values), [&] (const int i) {
                std::cout << i << " " << j << " " << local_mat.values(i,j) << std::endl;
            });
        });


        disp("kokkos_apply");
    }



        // each_read(P, [](const SizeType i, const SizeType j, const double value) {
        //     if(j != SizeType(value)) {
        //         std::cout << i << " " << j << " " << value << std::endl;
        //     }

        //    // utopia_test_assert(j == SizeType(value));

        // });





  
    void run_kokkos_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("KokkosTest");
        UTOPIA_RUN_TEST(kokkos_max);
        UTOPIA_RUN_TEST(kokkos_min);
        UTOPIA_RUN_TEST(kokkos_min_reduction);
        UTOPIA_RUN_TEST(kokkos_max_reduction);
        UTOPIA_RUN_TEST(kokkos_sum_reduction);
        UTOPIA_RUN_TEST(kokkos_write);
        UTOPIA_RUN_TEST(kokkos_read);
        UTOPIA_RUN_TEST(kokkos_apply)
        UTOPIA_UNIT_TEST_END("KokkosTest");
    }
}

#else //WITH_TRILINOS

namespace utopia
{
    void run_kokkos_test() {}
}

#endif //WITH_TRILINOS

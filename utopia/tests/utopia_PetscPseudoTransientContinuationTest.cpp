#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "utopia_assemble_laplacian_1D.hpp"
#include "test_problems/utopia_TestProblems.hpp"

namespace utopia
{


#ifdef WITH_PETSC
    class PseudoTransientContinuationTest {
    public:

        void run()
        {
            UTOPIA_RUN_TEST(PTC_small_test);
            UTOPIA_RUN_TEST(ASTRUM_small_test);

            // UTOPIA_RUN_TEST(affine_similarity_stiff_test);
        }


        void PTC_small_test()
        {
        	if(mpi_world_size() >1)
        		return;

        	SmallSingularExample<DMatrixd, DVectord> fun;
        	DVectord x_exact 	= values(2, 1.0);
        	DVectord x   		= values(2, 1.0);

        	{
        		Write<DVectord> r1(x_exact, LOCAL);
        		Write<DVectord> r2(x, LOCAL);

        		x.set(1, 0.0);
        		x_exact.set(0, 0.0);
        	}

            auto linear_solver = std::make_shared<GMRES<DMatrixd, DVectord> >(); 
            linear_solver->atol(1e-14); 
            linear_solver->max_it(10000);

            PseudoContinuation<DMatrixd, DVectord> solver(linear_solver); 
            solver.reset_mass_matrix(false); 

        	solver.verbose(verbose_);
        	solver.atol(1e-9);
        	solver.solve(fun, x);
        	utopia_test_assert(approxeq(x, x_exact, 1e-6));
        }

        void ASTRUM_small_test()
        {
            if(mpi_world_size() >1)
                return;

            SmallSingularExample<DMatrixd, DVectord> fun;
            DVectord x_exact    = values(2, 1.0);
            DVectord x          = values(2, 1.0);

            {
                Write<DVectord> r1(x_exact, LOCAL);
                Write<DVectord> r2(x, LOCAL);

                x.set(1, 0.0);
                x_exact.set(0, 0.0);
            }

            auto linear_solver = std::make_shared<GMRES<DMatrixd, DVectord> >(); 
            linear_solver->atol(1e-14); 
            linear_solver->max_it(10000);

            ASTRUM<DMatrixd, DVectord> solver(linear_solver); 
            solver.tau_init(2); 
            solver.scaling(false); 

            // solver.reset_mass_matrix(true); 

            solver.verbose(verbose_);
            solver.atol(1e-9);
            solver.solve(fun, x);
            utopia_test_assert(approxeq(x, x_exact, 1e-6));
        }        


        // void affine_similarity_stiff_test()
        // {
        // 	if(mpi_world_size() >1)
        // 		return;


        // 	const SizeType n = 100;

        // 	MildStiffExample<DMatrixd, DVectord> fun(n);
        // 	DVectord x, g;
        // 	fun.get_initial_guess(x);

        // 	auto linear_solver = std::make_shared<GMRES<DMatrixd, DVectord>>();
        // 	linear_solver->atol(1e-14);
        // 	linear_solver->max_it(10000);

        // 	AffineSimilarity<DMatrixd, DVectord> solver(linear_solver);

        // 	DMatrixd I = identity(n,n);
        // 	solver.set_mass_matrix(I);
        // 	solver.set_scaling_matrix(I);
        // 	solver.verbose(false);
        // 	solver.atol(1e-9);
        // 	solver.stol(1e-14);
        // 	solver.max_it(500);
        // 	solver.verbosity_level(VERBOSITY_LEVEL_NORMAL);
        // 	solver.solve(fun, x);
        // }

        PseudoTransientContinuationTest()
        : _n(100), verbose_(true) { }

    private:
        int _n;
        bool verbose_; 
    };

#endif //WITH_PETSC


    void runPetscPseudoTransientContinuationTest()
    {
        UTOPIA_UNIT_TEST_BEGIN("runPetscPseudoTransientContinuationTest");
        #ifdef WITH_PETSC
                PseudoTransientContinuationTest().run();
        #endif
        UTOPIA_UNIT_TEST_END("runPetscPseudoTransientContinuationTest");
    }
}

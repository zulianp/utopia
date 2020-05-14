#include "test_problems/utopia_TestProblems.hpp"
#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

namespace utopia
{


#ifdef WITH_PETSC
    class PseudoTransientContinuationTest {
    public:
        using Traits   = utopia::Traits<PetscVector>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

        void run()
        {
            // UTOPIA_RUN_TEST(PTC_small_test);
            // UTOPIA_RUN_TEST(ASTRUM_small_test);
            // UTOPIA_RUN_TEST(CSR_test_paper);
            // UTOPIA_RUN_TEST(affine_similarity_stiff_test);
        }


        void PTC_small_test()
        {
            // should use serial communicator
        	// if(mpi_world_size() >1)
        	// 	return;

        	SmallSingularExample<PetscMatrix, PetscVector> fun;
        	PetscVector x_exact(serial_layout(2), 1.0);
        	PetscVector x(serial_layout(2), 1.0);

        	{
        		Write<PetscVector> r1(x_exact, LOCAL);
        		Write<PetscVector> r2(x, LOCAL);

        		x.set(1, 0.0);
        		x_exact.set(0, 0.0);
        	}

            auto linear_solver = std::make_shared<GMRES<PetscMatrix, PetscVector> >();
            linear_solver->atol(1e-14);
            linear_solver->max_it(10000);

            PseudoContinuation<PetscMatrix, PetscVector> solver(linear_solver);
            solver.reset_mass_matrix(false);

        	solver.verbose(verbose_);
        	solver.atol(1e-9);
        	solver.solve(fun, x);
        	utopia_test_assert(approxeq(x, x_exact, 1e-6));
        }

        void ASTRUM_small_test()
        {
            // should use serial communicator
            // if(mpi_world_size() >1)
            //     return;

            SmallSingularExample<PetscMatrix, PetscVector> fun;
            PetscVector x_exact(serial_layout(2), 1.0);
            PetscVector x(serial_layout(2), 1.0);

            {
                Write<PetscVector> r1(x_exact, LOCAL);
                Write<PetscVector> r2(x, LOCAL);

                x.set(1, 0.0);
                x_exact.set(0, 0.0);
            }

            auto linear_solver = std::make_shared<GMRES<PetscMatrix, PetscVector> >();
            linear_solver->atol(1e-14);
            linear_solver->max_it(10000);

            ASTRUM<PetscMatrix, PetscVector> solver(linear_solver);
            solver.tau_init(2);
            solver.scaling(false);

            // solver.reset_mass_matrix(true);

            solver.verbose(verbose_);
            solver.atol(1e-9);
            solver.solve(fun, x);
            utopia_test_assert(approxeq(x, x_exact, 1e-6));
        }


        void CSR_test_paper()
        {
            // should use serial communicator
        	// if(mpi_world_size() >1)
        	// 	return;


        	ContinuousStirredReactor<PetscMatrix, PetscVector> fun;
        	PetscVector x0, x;
        	fun.get_initial_guess(x0, 0.0);

        	auto linear_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>(MATSOLVERPETSC, PCLU);
            // ASTRUM<PetscMatrix, PetscVector> solver(linear_solver);
            PseudoContinuation<PetscMatrix, PetscVector> solver(linear_solver);
            // solver.scaling(false);
            solver.reset_mass_matrix(true);
            solver.verbose(verbose_);
            solver.atol(1e-8);
            solver.max_it(5000000);

            // std::vector<double> vec_tau = {1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6};
            // std::vector<double> vec_tau = {1, 1e2, 1e4, 1e6};
            // std::vector<double> vec_tau = {1e-4, 1e-2};
            // std::vector<double> vec_tau = {1e-8};

            std::vector<double> vec_tau = {1e-4};

            for (double i : vec_tau) {
                x = x0;
                solver.tau_init(i);
                std::cout << "---- Solve with tau: " << i << " \n \n";
                solver.solve(fun, x);
            }

            // utopia_test_assert(approxeq(x, fun.exact_sol(), 1e-4));
        }

        PseudoTransientContinuationTest() = default;

    private:
        int _n{100};
        bool verbose_{true};
    };

#endif //WITH_PETSC


    static void pseudo_transient_continuation()
    {
#ifdef WITH_PETSC
        PseudoTransientContinuationTest().run();
#endif
    }

    UTOPIA_REGISTER_TEST_FUNCTION(pseudo_transient_continuation);
}  // namespace utopia

#ifndef UTOPIA_BENCHMARK_ALGORITHMS_HPP
#define UTOPIA_BENCHMARK_ALGORITHMS_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_assemble_laplacian_1D.hpp"
#include "utopia_TestProblems.hpp"
#include "utopia_BiCGStab.hpp"

#include <string>
#include <cassert>

namespace utopia {

    template<class Matrix, class Vector>
    class BenchmarkAlgorithms : public Benchmark {
    public:
        using Traits   = utopia::Traits<Vector>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

        BenchmarkAlgorithms(const Comm &comm = Comm()) : comm_(comm) {}


        virtual std::string name() override
        {
            return "Algorithms";
        }

        void initialize() override
        {
            static const bool is_sparse = utopia::is_sparse<Matrix>::value;
            static const bool verbose = false;

            const SizeType base_n = is_sparse? 1000 : 10;
            const SizeType n_instances = 5;


            for(SizeType i = 0; i < n_instances; ++i) {
                const SizeType n = base_n * (i + 1);

                //Conjugate gradient method
                this->register_experiment(
                    "cg_" + std::to_string(i),
                    [n]() {
                        ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
                        cg.verbose(verbose);
                        cg.max_it(n * mpi_world_size());
                        cg.set_preconditioner(std::make_shared< InvDiagPreconditioner<Matrix, Vector> >());
                        run_linear_solver(n, cg);
                    }
                );

                // this->register_experiment(
                //     "gs_" + std::to_string(i),
                //     [n]() {
                //         GaussSeidel<Matrix, Vector, HOMEMADE> gs;
                //         gs.verbose(true);
                //         gs.max_it(100*n * mpi_world_size());
                //         gs.check_convergence_each(50);
                //         run_linear_solver(n, gs);
                //     }
                // );

                // this->register_experiment(
                //     "l1_gs_" + std::to_string(i),
                //     [n]() {
                //         GaussSeidel<Matrix, Vector, HOMEMADE> gs;
                //         gs.l1(true);
                //         gs.verbose(true);
                //         gs.max_it(100*n * mpi_world_size());
                //         gs.check_convergence_each(50);
                //         run_linear_solver(n, gs);
                //     }
                // );

                this->register_experiment(
                    "bicgstab_" + std::to_string(i),
                    [n]() {
                        BiCGStab<Matrix, Vector, HOMEMADE> cg;
                        cg.verbose(verbose);
                        cg.max_it(n * mpi_world_size());
                        run_linear_solver(n, cg);
                    }
                );

                this->register_experiment(
                    "newton_cg_" + std::to_string(i),
                    [i]() {
                        Rastrigin<Matrix, Vector> fun(10 * (i+1));
                        Vector x = local_values(10 * (i+1), 1.);

                        ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
                        cg.max_it(size(x).get(0));

                        auto backtracking = std::make_shared<utopia::Backtracking<Vector> >();

                        Newton<Matrix, Vector, HOMEMADE> newton(make_ref(cg));
                        newton.set_line_search_strategy(backtracking);

                        double mag_x0 = -1;
                        fun.value(x, mag_x0);

                        newton.solve(fun, x);

                        double mag_x = -1.;
                        fun.value(x, mag_x);
                        utopia_test_assert(mag_x <= mag_x0);
                    }
                );

                this->register_experiment(
                    "trust_region_" + std::to_string(i),
                    [i]() {
                        Rastrigin<Matrix, Vector> fun(10 * (i+1));
                        Vector x = local_values(10 * (i+1), 1.);

                        auto st_cg = std::make_shared<SteihaugToint<Matrix, Vector> >();

                        TrustRegion<Matrix, Vector> trust_region(st_cg);
                        trust_region.verbose(false);

                        double mag_x0 = -1;
                        fun.value(x, mag_x0);

                        trust_region.solve(fun, x);

                        double mag_x = -1.;
                        fun.value(x, mag_x);
                        utopia_test_assert(mag_x <= mag_x0);
                    }
                );

                this->register_experiment(
                    "projected_gradient_" + std::to_string(i),
                    [i]() {
                        ProjectedGradient<Matrix, Vector, HOMEMADE> pg;
                        run_qp_solver((base_n/2) * (i + 1), pg);
                    }
                );

                this->register_experiment(
                    "projected_conjugate_gradient_" + std::to_string(i),
                    [i]() {
                        ProjectedConjugateGradient<Matrix, Vector, HOMEMADE> pg;
                        run_qp_solver((base_n/2) * (i + 1), pg);
                    }
                );

                this->register_experiment(
                    "projected_gauss_seidel_" + std::to_string(i),
                    [i]() {
                        ProjectedGaussSeidel<Matrix, Vector, HOMEMADE> pg;
                        run_qp_solver((base_n/2) * (i + 1), pg);
                    }
                );

                this->register_experiment(
                    "projected_l1_gauss_seidel_" + std::to_string(i),
                    [i]() {
                        ProjectedGaussSeidel<Matrix, Vector, HOMEMADE> pg;
                        // pg.l1(true);
                        run_qp_solver((base_n/2) * (i + 1), pg);
                    }
                );

                //FIXME
                // this->register_experiment("multigrid_" + std::to_string(i), [n]() {

                // 	auto smoother      = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
                // 	auto coarse_solver = std::make_shared<BiCGStab<Matrix, Vector, HOMEMADE>>();

                // 	// smoother->set_preconditioner(std::make_shared< InvDiagPreconditioner<Matrix, Vector> >());
                // 	coarse_solver->set_preconditioner(std::make_shared< InvDiagPreconditioner<Matrix, Vector> >());
                // 	coarse_solver->verbose(true);
                // 	coarse_solver->max_it(2000);
                // 	// coarse_solver->reset_initial_guess(true);
                // 	coarse_solver->atol(1e-16);

                // 	Multigrid<Matrix, Vector, HOMEMADE> multigrid(
                // 	                                    smoother,
                // 	                                    coarse_solver
                // 	                                    );

                // 	run_multigrid(n, multigrid);
                // });

            }
        }

    private:

        template<class QPSolver>
        static void run_qp_solver(const Comm &comm, const SizeType n, QPSolver &qp_solver)
        {
            auto vl = layout(comm, n, n * comm.size());
            auto ml = layout(comm, n, n, n * comm.size(), n * comm.size());

            Matrix m; m.sparse(ml, 3, 2);
            assemble_laplacian_1D(m);

            auto N = size(m).get(0);

            {
                Range r = row_range(m);
                Write<Matrix> w(m);
                if(r.inside(0)) {
                    m.set(0, 0, 1.);
                    m.set(0, 1, 0);
                }

                if(r.inside(N)) {
                    m.set(N-1, N-1, 1.);
                    m.set(N-1, N-2, 0);
                }
            }

            Vector rhs(vl, 1.);
            {
                //Creating test vector (alternative way see [assemble vector alternative], which might be easier for beginners)
                Range r = range(rhs);
                Write<Vector> w(rhs);

                if(r.inside(0)) {
                    rhs.set(0, 0);
                }

                if(r.inside(N-1)) {
                    rhs.set(N-1, 0.);
                }
            }

            Vector upper_bound(vl, 100.0);
            Vector solution(vl, 0.0);

            qp_solver.max_it(N*2);
            // qp_solver.verbose(true);
            qp_solver.set_box_constraints(make_upper_bound_constraints(make_ref(upper_bound)));

            bool ok = qp_solver.solve(m, rhs, solution);
            utopia_test_assert(ok);
        }

        static void run_linear_solver(const Comm &comm, const SizeType n, LinearSolver<Matrix, Vector> &solver)
        {
            auto vl = layout(comm, n, n * comm.size());
            auto ml = layout(comm, n, n, n * comm.size(), n * comm.size());

            Matrix A; A.sparse(ml, 3, 2);
            Vector b(vl, 1.);
            Vector x(vl, 0.);

            assemble_laplacian_1D(A, true);

            auto N = size(A).get(0);
            {
                Range r = row_range(A);
                Write<Vector> w_b(b);

                if(r.inside(0)) {
                    b.set(0, 0.);
                }

                if(r.inside(N-1)) {
                    b.set(N-1, 0.);
                }
            }

            solver.solve(A, b, x);
            Vector Ax = A * x;

            const Scalar r_norm = norm2(b - Ax);

            utopia_test_assert(approxeq(Ax, b, 1e-6));
        }

        template<class MultigridSolver>
        static void run_multigrid(const SizeType n, MultigridSolver &multigrid)
        {
            using TransferT       = utopia::Transfer<Matrix, Vector>;
            using IPTransferT     = utopia::IPTransfer<Matrix, Vector>;
            using MatrixTransferT = utopia::MatrixTransfer<Matrix, Vector>;

            const static bool verbose   = false;
            const static bool use_masks = true;
            const SizeType n_levels = 5;

            MultiLevelTestProblem1D<Matrix, Vector, Poisson1D<Matrix, Vector> > ml_problem(n/pow(2, n_levels-1), n_levels, !use_masks);

            multigrid.max_it(50);
            multigrid.atol(1e-13);
            multigrid.stol(1e-13);
            multigrid.rtol(1e-9);
            multigrid.pre_smoothing_steps(3);
            multigrid.post_smoothing_steps(3);
            multigrid.set_fix_semidefinite_operators(true);
            multigrid.must_generate_masks(use_masks);;
            multigrid.verbose(verbose);

            auto transfers = ml_problem.get_transfer();
            auto funs = ml_problem.get_functions();

            multigrid.set_transfer_operators(transfers);

            Vector x, g;
            funs[n_levels - 1]->get_eq_constrains_values(x);

            Matrix A;
            funs[n_levels - 1]->hessian(x, A);
            funs[n_levels - 1]->gradient(x, g);
            multigrid.update(std::make_shared<Matrix>(A));

            if(verbose) {
                multigrid.describe();
            }

            // x = *ml_problem.rhs;
            multigrid.apply(g, x);

            double diff0 = norm2(A * x);
            double diff  = norm2(g - A * x);
            double rel_diff = diff/diff0;

            utopia_test_assert(rel_diff < 1e-8);
        }

    private:
        Comm comm_;
    };
}

#endif //UTOPIA_BENCHMARK_ALGORITHMS_HPP

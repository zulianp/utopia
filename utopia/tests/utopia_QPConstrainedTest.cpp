#include <cassert>
#include <string>
#include "utopia.hpp"
#include "utopia_Base.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_ConstrainedBenchmark.hpp"
#include "utopia_LargeScaleIncludes.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Testing.hpp"

namespace utopia {

    // FIXME Move this to the benchmark target (test should not have output unless they fail)
    template <class Matrix, class Vector>
    class QPConstrainedBenchmark final : public Benchmark {
    public:
        DEF_UTOPIA_SCALAR(Vector);

        std::string name() override { return "QPConstrainedBenchmark benchmark."; }

        explicit QPConstrainedBenchmark(const SizeType &n, const bool verbose) : n_(n), verbose_(verbose) {
            test_functions_.resize(7);
            test_functions_[0] = std::make_shared<Poisson1D<Matrix, Vector> >(n_ * mpi_world_size(), 1);
            test_functions_[1] = std::make_shared<Poisson1D<Matrix, Vector> >(n_ * mpi_world_size(), 2);
            test_functions_[2] = std::make_shared<Poisson1D<Matrix, Vector> >(n_ * mpi_world_size(), 3);
            test_functions_[3] = std::make_shared<Poisson1D<Matrix, Vector> >(n_ * mpi_world_size(), 4);

            // // works only with  petsc
            test_functions_[4] = std::make_shared<Poisson2D<Matrix, Vector> >(n_ * mpi_world_size(), 1);
            test_functions_[5] = std::make_shared<Poisson2D<Matrix, Vector> >(n_ * mpi_world_size(), 2);
            test_functions_[6] = std::make_shared<Membrane2D<Matrix, Vector> >(n_ * mpi_world_size());
        }

        ~QPConstrainedBenchmark() override { test_functions_.clear(); }

        void initialize() override {
            this->register_experiment("MPGRP_Test", [this]() {
                MPGRP<Matrix, Vector> solver;

                solver.verbose(true);
                run_test(this->test_functions_, solver, "MPGRP_Test", this->verbose_);
            });

            this->register_experiment("ProjectedGradient_Test", [this]() {
                ProjectedGradient<Matrix, Vector> solver;

                solver.verbose(true);
                run_test(this->test_functions_, solver, "ProjectedGradient_Test", this->verbose_);
            });

            this->register_experiment("ProjectedGaussSeidel_Test", [this]() {
                ProjectedGaussSeidel<Matrix, Vector> solver;

                solver.verbose(true);
                run_test(this->test_functions_, solver, "ProjectedGaussSeidel_Test", this->verbose_);
            });

            // has problems with some tests - to be debugged
            // this->register_experiment("ProjectedConjugateGradient_Test",
            // 	[this]() {
            // 		ProjectedConjugateGradient<Matrix, Vector> solver;

            //            solver.verbose(true);
            //            run_test(this->test_functions_, solver, "ProjectedConjugateGradient_Test", this->verbose_);
            // 	}
            // );

            // works only for petsc
#ifdef UTOPIA_WITH_PETSC
            this->register_experiment("ProjectedTao_Test", [this]() {
                auto lin_solver = std::make_shared<GMRES<Matrix, Vector> >();
                TaoQPSolver<Matrix, Vector> solver(lin_solver);
                solver.tao_type("gpcg");
                solver.verbose(true);

                run_test(this->test_functions_, solver, "ProjectedTao_Test", this->verbose_);
            });
#endif  // UTOPIA_WITH_PETSC
        }

    private:
        template <class Fun, class NonlinearSolver>
        static void run_test(std::vector<std::shared_ptr<Fun> > &test_functions,
                             NonlinearSolver &solver,
                             const std::string &solv_name,
                             const bool &exp_verbose) {
            if (exp_verbose && mpi_world_rank() == 0) {
                utopia::out() << "--------------------------------------------------------- \n";
                utopia::out() << "				" << solv_name << "				\n";
                utopia::out() << "--------------------------------------------------------- \n";
            }

            InputParameters in;
            in.set("atol", 1e-6);
            in.set("rtol", 1e-11);
            in.set("stol", 1e-14);
            in.set("stol", 1e-14);
            in.set("delta_min", 1e-13);
            in.set("max-it", 10000);
            in.set("verbose", false);
            solver.read(in);

            for (size_t i = 0; i < test_functions.size(); i++) {
                Vector x_init = test_functions[i]->initial_guess();
                auto lo = layout(x_init);

                Vector g;
                Matrix H;
                test_functions[i]->gradient(x_init, g);
                g *= -1.0;
                test_functions[i]->hessian(x_init, H);

                auto box = test_functions[i]->box_constraints();
                box.fill_empty_bounds(lo);

                auto lb = std::make_shared<Vector>();
                lb->zeros(lo);

                auto ub = std::make_shared<Vector>();
                ub->zeros(lo);

                auto correction_constraints = make_box_constaints(lb, ub);

                *correction_constraints.upper_bound() = *box.upper_bound() - x_init;
                *correction_constraints.lower_bound() = *box.lower_bound() - x_init;

                Vector s;
                s.zeros(lo);
                solver.set_box_constraints(correction_constraints);
                solver.solve(H, g, s);

                // taking correction
                x_init += s;

                bool feas_flg = test_functions[i]->is_feasible(x_init);
                utopia_test_assert(feas_flg);

                // auto sol_status = solver.solution_status();

                // Membrane2D<Matrix, Vector> * fun_poisson2D = dynamic_cast<Membrane2D<Matrix, Vector>
                // *>(test_functions[i].get()); fun_poisson2D->output_to_VTK(x_init, "Membrane2D.vtk");

                // // disp(x_init);

                // Scalar energy;
                // test_functions[i]->value(x_init, energy);
                // std::cout<<"energy: "<< energy << "  \n";

                // if(exp_verbose && mpi_world_rank()==0)
                // {
                // 	std::cout<< i <<std::setw(5-std::to_string(i).size()) <<" : "<< test_functions[i]->name() <<"_"
                // << dim <<  std::right <<  std::setw(60-std::to_string(dim).size() - test_functions[i]->name().size())
                // << std::right << "its:  " << num_its << std::setw(5-std::to_string(num_its).size())<<  "  \n";

                // 	const auto conv_reason = sol_status.reason;
                // 	if(conv_reason< 0)
                // 	{
                // 		sol_status.describe(std::cout);
                // 	}

                // }
            }
        }

    private:
        std::vector<std::shared_ptr<ConstrainedExtendedTestFunction<Matrix, Vector> > > test_functions_;

        SizeType n_;
        bool verbose_;
    };

    static void qp_constrained() {
#ifdef UTOPIA_WITH_PETSC
        int verbosity_level = 1;
        const int n_global = 20;
        bool alg_verbose = false;

        if (Utopia::instance().verbose()) {
            verbosity_level = 2;
        }

        QPConstrainedBenchmark<PetscMatrix, PetscVector> bench_petsc(n_global, alg_verbose);
        bench_petsc.set_verbosity_level(verbosity_level);
        bench_petsc.run();
#endif  // UTOPIA_WITH_PETSC

        // #ifdef WITH_TRILINOS
        //     QPConstrainedBenchmark<TpetraMatrixd, TpetraVectord> bench_tril(n_global, alg_verbose);
        // 	bench_tril.set_verbosity_level(verbosity_level);
        // 	bench_tril.run();
        // #endif //WITH_TRILINOS

        // #ifdef WITH_BLAS
        // 	QPConstrainedBenchmark<BlasMatrixd, BlasVectord> bench_blas(n_global, alg_verbose);
        // 	bench_blas.set_verbosity_level(verbosity_level);
        // 	bench_blas.run();
        // #endif //WITH_BLAS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(qp_constrained);
}  // namespace utopia

// #include "utopia_trilinos_Base.hpp"
#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
// #include "utopia_trilinos_DeviceView.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_For.hpp"
#include "utopia_Allocations.hpp"

namespace utopia
{

    template<typename Matrix, typename Vector>
    class HckTests
    {

    public:
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef UTOPIA_SCALAR(Vector) Scalar;

        HckTests(const SizeType & n, const SizeType & n_levels=2, const Scalar & lambda = 1.0, const bool verbose = false, const bool output_flg = false):
        n_(n),
        n_levels_(n_levels),
        lambda_(lambda),
        verbose_(verbose),
        output_vtk_(output_flg)
        {
            input_params_.set("atol", 1e-5);
            input_params_.set("rtol", 1e-10);
            input_params_.set("stol", 1e-10);
            input_params_.set("verbose", verbose_);
            input_params_.set("max-it", 20);

            // RMTR specific parameters
            input_params_.set("max_coarse_it", 2);
            input_params_.set("max_sucessful_coarse_it", 1);
            input_params_.set("max_QP_coarse_it", 1000);
            input_params_.set("pre_smoothing_steps", 5);
            input_params_.set("post_smoothing_steps", 5);
            input_params_.set("max_sucessful_smoothing_it", 1);
            input_params_.set("max_QP_smoothing_it", 10);
            input_params_.set("delta0", 1e3);
            input_params_.set("grad_smoothess_termination", 1e-8);
        }

        void run_petsc()
        {
            UTOPIA_RUN_TEST(newton_test);

            UTOPIA_RUN_TEST(STCG_test);
            UTOPIA_RUN_TEST(MPGRP_test);
            UTOPIA_RUN_TEST(ProjectedGS);


            UTOPIA_RUN_TEST(TR_unconstrained);
            UTOPIA_RUN_TEST(TR_constrained);

            UTOPIA_RUN_TEST(Poisson_test);

            UTOPIA_RUN_TEST(QuasiTR_unconstrained);
            UTOPIA_RUN_TEST(QuasiTR_constrained);


            UTOPIA_RUN_TEST(RMTR_unconstrained);
            UTOPIA_RUN_TEST(RMTR_l2_linear);

            UTOPIA_RUN_TEST(RMTR_inf_linear_unconstr);

        }

        void run_trilinos()
        {
            // UTOPIA_RUN_TEST(TR_tril_test);

            // UTOPIA_RUN_TEST(RMTR_l2_test);
            // UTOPIA_RUN_TEST(RMTR_inf_test);
            // UTOPIA_RUN_TEST(Quasi_RMTR_l2_test);

            // UTOPIA_RUN_TEST(Quasi_RMTR_inf_test);
            UTOPIA_RUN_TEST(transform_test);
            UTOPIA_RUN_TEST(e_mul_test);
            UTOPIA_RUN_TEST(e_div_test);
            UTOPIA_RUN_TEST(negate_alpha_test);
            UTOPIA_RUN_TEST(axpy_test);
            UTOPIA_RUN_TEST(quad_form_test);


           // UTOPIA_RUN_TEST(STCG_test);
            // UTOPIA_RUN_TEST(for_each_loop_test);
            // UTOPIA_RUN_TEST(parallel_each_write_test);

            // UTOPIA_RUN_TEST(quad_form_test);

            // UTOPIA_RUN_TEST(multi_reduce_test);
            UTOPIA_RUN_TEST(MPGRP_test);


            // UTOPIA_RUN_TEST(residual_test);


            // //THIS
            // UTOPIA_RUN_TEST(Quasi_RMTR_inf_test);
        }

        void e_mul_test()
        {
            Vector x = values(n_, 1.0);
            Vector y = values(n_, 2.0);
            Vector z = values(n_, 0.0);

            UTOPIA_NO_ALLOC_BEGIN("e_mul_test");
            z = e_mul(x, y);
            UTOPIA_NO_ALLOC_END();
        }

        void e_div_test()
        {
            Vector x = values(n_, 6.0);
            Vector z = values(n_, 3.0);

            UTOPIA_NO_ALLOC_BEGIN("e_div_test");
            z = x / z;
            
            Scalar sum_z = sum(z);
            utopia_test_assert(approxeq(sum_z, 2.0*n_));
            z = x / x;
            sum_z = sum(z);
            utopia_test_assert(approxeq(sum_z, 1.0*n_));

            z = z / x;
            sum_z = sum(z);
            utopia_test_assert(approxeq(sum_z, 1.0/6.0*n_));

            UTOPIA_NO_ALLOC_END();
        }

        void quad_form_test()
        {
            Vector x = values(n_, 1.0);
            Vector y = values(n_, 2.0);

            auto expr = 0.5 * dot(x, y) - 0.5 * dot(x, y);

            Scalar val = expr;

            utopia_test_assert(approxeq(val, 0.0));
        }

        void residual_test()
        {
            Vector x = values(n_, 1.0);
            Vector b = values(n_, 1.0);
            Matrix A = sparse(n_, n_, 3);
            Vector r = values(n_, 0.0);

            UTOPIA_NO_ALLOC_BEGIN("residual_test");
            r = x - b;
            UTOPIA_NO_ALLOC_END();
        }

        void axpy_test()
        {
            Vector x = values(n_, 1.0);
            Vector y = values(n_, 1.0);
            Vector p = values(n_, 2.0);
            
            UTOPIA_NO_ALLOC_BEGIN("axpy_test");
            x = x - 0.5 * p;
            y = x - 0.5 * p;
            UTOPIA_NO_ALLOC_END();
        }

        void transform_test()
        {
            Vector x = values(n_, 1.0);
            
            parallel_transform(
                x,
                UTOPIA_LAMBDA(const SizeType &i, const Scalar &v) -> Scalar {
                    return (i+1)*v;
            });
            
            Scalar expected = ((n_ + 1) * n_)/2.0;
            Scalar sum_x = sum(x);

            utopia_test_assert(approxeq(sum_x, expected, 1e-10));
        }

        void negate_alpha_test()
        {
            Vector x = values(n_, 1.0);
            Vector y = values(n_, 1.0);
            
            UTOPIA_NO_ALLOC_BEGIN("negate_alpha_test");
            y = -0.5 * x;
            UTOPIA_NO_ALLOC_END();
        }


        template<class QPSolverTemp>
        void QP_solve(QPSolverTemp &qp_solver) const
        {
            Matrix H_working;
            Vector g_working, x_working;

#ifdef WITH_PETSC
            Bratu2D<PetscMatrix, PetscVector> fun(n_);
            PetscVector x = fun.initial_guess();
            PetscMatrix H;
            PetscVector g;

            fun.hessian(x, H);
            fun.gradient(x, g);

            backend_convert_sparse(H, H_working);
            backend_convert(g, g_working);
            x_working =  0.0 * g_working;

            // monitor(0, H, "Hessian.m", "H");
            // monitor(0, g, "gradient.m", "g");

            if(dynamic_cast<QPSolver<Matrix, Vector> *>(qp_solver.get()))
            {
                QPSolver<Matrix, Vector> * qp_box = dynamic_cast<QPSolver<Matrix, Vector> *>(qp_solver.get());
                Vector lb = local_values(local_size(x_working).get(0), -9e9);
                Vector ub = local_values(local_size(x_working).get(0), 9e9);
                qp_box->set_box_constraints(make_box_constaints(make_ref(lb), make_ref(ub)));
                qp_box->solve(H_working, -1.0*g_working, x_working);
            }
            else
            {
                qp_solver->solve(H_working, -1.0*g_working, x_working);
            }

#endif

        }


        void STCG_test()
        {

            // auto QP_solver = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "sor", false);
            auto QP_solver = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
            QP_solver->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());
            // auto precond = std::make_shared<GaussSeidel<Matrix, Vector, HOMEMADE> >();
            // precond->verbose(verbose_);
            // precond->max_it(1);
            // precond->use_line_search(false);
            // QP_solver->set_preconditioner(precond);
            QP_solver->use_precond_direction(false);



            QP_solver->atol(1e-10);
            // QP_solver->max_it(n_*n_);
            QP_solver->max_it(10);
            QP_solver->verbose(verbose_);
            // QP_solver->norm_schedule(NormSchedule::NEVER); 
            QP_solver->norm_schedule(NormSchedule::EVERY_ITER); 
            QP_solver->current_radius(1e5);

            QP_solve(QP_solver);
        }

        void MPGRP_test()
        {
            auto QP_solver = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
            QP_solver->atol(1e-10);
            QP_solver->max_it(2);
            QP_solver->verbose(true);

            QP_solve(QP_solver);
        }

        void ProjectedGS()
        {
            auto QP_solver = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();
            QP_solver->atol(1e-10);
            QP_solver->max_it(n_*n_);
            QP_solver->verbose(verbose_);

            QP_solve(QP_solver);
        }

        void Poisson_test()
        {
//FIXME
#ifdef WITH_PETSC
            Poisson3D<Matrix, Vector> fun(10);

            Vector b;
            Matrix H;
            Vector x = fun.initial_guess();

            if(verbose_)
                fun.describe();

#ifdef WITH_PETSC
            auto subproblem = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", false);
#else
            // auto subproblem = std::make_shared<utopia::Lanczos<Matrix, Vector> >();
            auto subproblem = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
#endif
            // subproblem->pc_type("bjacobi");
            subproblem->atol(1e-14);
            subproblem->max_it(1000);

            TrustRegion<Matrix, Vector> tr_solver(subproblem);
            tr_solver.read(input_params_);
            tr_solver.delta0(0.01);
            tr_solver.atol(1e-12);
            tr_solver.solve(fun, x);


            PetscMultilevelTestProblem<Matrix, Vector, Poisson3D<Matrix, Vector> > multilevel_problem(3, n_levels_, n_);

            auto fun1 = multilevel_problem.level_functions_[n_levels_-1];
            Poisson3D<Matrix, Vector> * fun_Laplace = dynamic_cast<Poisson3D<Matrix, Vector> *>(fun1.get());


            x = fun_Laplace->initial_guess();
            fun_Laplace->gradient(x, b);
            fun_Laplace->hessian(x, H);


            auto direct_solver = std::make_shared<utopia::RedundantLinearSolver<Matrix, Vector> >();
            direct_solver->number_of_parallel_solves(mpi_world_size());

            auto smoother = std::make_shared<GaussSeidel<PetscMatrix, PetscVector>>();

            Multigrid<Matrix, Vector> multigrid(smoother, direct_solver);
            multigrid.set_transfer_operators(multilevel_problem.transfers_);
            multigrid.update(make_ref(H));
            multigrid.read(input_params_);
            multigrid.apply(b, x);


            if(output_vtk_)
                fun.output_to_VTK(x);
#endif //WITH_PETSC
        }


        void TR_unconstrained()
        {
            Bratu2D<Matrix, Vector> fun(100, 5.0);
            Vector x = fun.initial_guess();

            if(verbose_)
                fun.describe();

#ifdef WITH_PETSC
            auto subproblem = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", false);
#else
            auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector>>();
#endif

            TrustRegion<Matrix, Vector> tr_solver(subproblem);
            tr_solver.read(input_params_);
            tr_solver.verbose(verbose_);
            tr_solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x);
        }


        void newton_test()
        {
            Bratu2D<Matrix, Vector> fun(n_, 5.0);
            Vector x = fun.initial_guess();

            if(verbose_)
                fun.describe();

            auto lsolver = std::make_shared<GMRES<Matrix, Vector> >();
            // lsolver->pc_type("bjacobi");

            Newton<Matrix, Vector> solver(lsolver);
            solver.read(input_params_);
            solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x);
        }

        void QuasiTR_unconstrained()
        {
            Bratu2D<Matrix, Vector> fun(n_);
            Vector x = fun.initial_guess();
            SizeType memory_size = 5;

            auto subproblem = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
            subproblem->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
            subproblem->atol(1e-14);
            subproblem->max_it(100000);

            auto hess_approx   = std::make_shared<LBFGS<Vector> >(memory_size);
            hess_approx->theta_min(1.0);
            hess_approx->damping_tech(POWEL);
            hess_approx->scaling_tech(ADAPTIVE);

            QuasiTrustRegion<Vector> tr_solver(hess_approx, subproblem);
            tr_solver.read(input_params_);
            tr_solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x, "QuasiTRUnstrained.vtk");
        }


        void TR_constrained()
        {
            Bratu2D<Matrix, Vector> fun(100, 5.0);
            Vector x = fun.initial_guess();

            if(verbose_)
                fun.describe();

            // auto qp_solver = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
            auto qp_solver = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();
            qp_solver->atol(1e-10);
            qp_solver->max_it(n_*n_);
            qp_solver->use_line_search(false);
            qp_solver->verbose(false);


            TrustRegionVariableBound<Matrix, Vector> tr_solver(qp_solver);

            Vector ub, lb;
            fun.upper_bound(ub);
            fun.lower_bound(lb);
            auto box = make_box_constaints(make_ref(lb), make_ref(ub));

            tr_solver.set_box_constraints(box);
            tr_solver.read(input_params_);
            tr_solver.verbose(verbose_);
            tr_solver.delta0(1e-2);
            tr_solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x);
        }


        void QuasiTR_constrained()
        {
            Bratu2D<Matrix, Vector> fun(n_);
            Vector x = fun.initial_guess();
            SizeType memory_size = 10;

            auto qp_solver = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
            qp_solver->atol(1e-10);
            qp_solver->max_it(n_*n_);

            auto hess_approx   = std::make_shared<LBFGS<Vector> >(memory_size);
            // hess_approx->theta_min(1.0);
            // hess_approx->damping_tech(POWEL);
            // hess_approx->scaling_tech(ADAPTIVE);

            Vector ub, lb;
            fun.upper_bound(ub);
            fun.lower_bound(lb);
            auto box = make_box_constaints(make_ref(lb), make_ref(ub));


            QuasiTrustRegionVariableBound<Vector>  tr_solver(hess_approx, qp_solver);
            tr_solver.set_box_constraints(box);
            tr_solver.read(input_params_);
            tr_solver.delta0(0.01);
            tr_solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x, "QuasiTRConstrained.vtk");
        }


        void RMTR_unconstrained()
        {
//FIXME
#ifdef WITH_PETSC
            PetscMultilevelTestProblem<Matrix, Vector, Bratu2D<Matrix, Vector> > multilevel_problem(2, n_levels_, n_);
            auto fun = multilevel_problem.level_functions_[n_levels_-1];
            Bratu2D<Matrix, Vector> * fun_Bratu2D = dynamic_cast<Bratu2D<Matrix, Vector> *>(fun.get());
            Vector x = fun_Bratu2D->initial_guess();

            if(verbose_)
                fun_Bratu2D->describe();

#ifdef WITH_PETSC
            auto tr_strategy_coarse = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", true);
#else
            auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
#endif
            auto tr_strategy_fine = std::make_shared<utopia::Lanczos<Matrix, Vector> >("sor");


            // auto rmtr = std::make_shared<RMTR<Matrix, Vector, SECOND_ORDER> >(n_levels_);
            auto rmtr = std::make_shared<RMTR<Matrix, Vector, SECOND_ORDER> >(n_levels_);

            // Set TR-QP strategies
            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);

            // Transfers and objective functions
            rmtr->set_transfer_operators(multilevel_problem.transfers_);
            rmtr->set_functions( multilevel_problem.level_functions_);


            rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);
            rmtr->read(input_params_);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);

            // Solve
            rmtr->solve(x);


            if(output_vtk_)
                fun_Bratu2D->output_to_VTK(x, "RMTR_output.vtk");
#endif //WITH_PETSC
        }

        void RMTR_l2_linear()
        {
#ifdef WITH_PETSC
            PetscMultilevelTestProblem<Matrix, Vector, Poisson3D<Matrix, Vector> > multilevel_problem(3, n_levels_, n_);

            auto fun = multilevel_problem.level_functions_[n_levels_-1];
            Poisson3D<Matrix, Vector> * fun_Poisson3D = dynamic_cast<Poisson3D<Matrix, Vector> *>(fun.get());
            Vector x = fun_Poisson3D->initial_guess();

            if(verbose_)
                fun_Poisson3D->describe();


            auto tr_strategy_coarse = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", true);
            auto tr_strategy_fine = std::make_shared<utopia::Lanczos<Matrix, Vector> >("sor");


            // auto rmtr = std::make_shared<RMTR<Matrix, Vector, GALERKIN> >(n_levels_);
            // auto rmtr = std::make_shared<RMTR<Matrix, Vector, SECOND_ORDER> >(n_levels_);
            auto rmtr = std::make_shared<RMTR<Matrix, Vector, FIRST_ORDER> >(n_levels_);

            // Set TR-QP strategies
            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);

            // Transfers and objective functions
            rmtr->set_transfer_operators(multilevel_problem.transfers_);
            rmtr->set_functions( multilevel_problem.level_functions_);


            rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
            rmtr->read(input_params_);

            // Solve
            rmtr->solve(x);

            if(output_vtk_)
                fun_Poisson3D->output_to_VTK(x, "RMTR__linear_output.vtk");
#endif //WITH_PETSC
        }


        void RMTR_inf_linear_unconstr()
        {
            #ifdef WITH_PETSC
                PetscMultilevelTestProblem<Matrix, Vector, Poisson3D<Matrix, Vector> > multilevel_problem(3, n_levels_, n_);

                auto fun = multilevel_problem.level_functions_[n_levels_-1];
                Poisson3D<Matrix, Vector> * fun_Poisson3D = dynamic_cast<Poisson3D<Matrix, Vector> *>(fun.get());
                Vector x = fun_Poisson3D->initial_guess();

                if(verbose_)
                    fun_Poisson3D->describe();

                // ---------------------- TODO:: investigate why we get negative alpha --------------
                // auto tr_strategy_coarse = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();
                // auto tr_strategy_fine = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();

                auto tr_strategy_coarse = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
                auto tr_strategy_fine   = std::make_shared<utopia::MPGRP<Matrix, Vector> >();

                auto rmtr = std::make_shared<RMTR_inf<Matrix, Vector, FIRST_ORDER> >(n_levels_);

                // Set TR-QP strategies
                rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
                rmtr->set_fine_tr_strategy(tr_strategy_fine);

                // Transfers and objective functions
                rmtr->set_transfer_operators(multilevel_problem.transfers_);
                rmtr->set_functions( multilevel_problem.level_functions_);


                rmtr->read(input_params_);
                rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);
                rmtr->verbose(verbose_);
                // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
                rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);


                // Solve
                rmtr->solve(x);

                if(output_vtk_)
                    fun_Poisson3D->output_to_VTK(x, "RMTR__linear_output.vtk");
            #endif //WITH_PETSC
        }


        void TR_tril_test()
        {
            #ifdef WITH_PETSC
                Poisson3D<PetscMatrix, PetscVector> fun(50);
                PetscVector x = fun.initial_guess();

                PetscMatrix H;
                PetscVector g, x_eq, x_bc_marker, empty_rhs;

                fun.get_A_rhs(H, g);
                fun.get_eq_constrains_values(x_eq);
                fun.get_eq_constrains_flg(x_bc_marker);

                empty_rhs = 0.0*x;


                PetscVector lb = local_values(local_size(x), -9e9);
                PetscVector ub = local_values(local_size(x), 9e9);

                QuadraticExtendedFunction<PetscMatrix, PetscVector> fun_QP(H, g, x_eq, x_bc_marker, empty_rhs);

                #ifdef WITH_TRILINOS
                    Matrix H_tril;
                    Vector g_tril;
                    Vector x_tril, x_eq_tril, x_bc_marker_tril, empty_rhs_tril;
                    backend_convert_sparse(H, H_tril);
                    backend_convert(g, g_tril);
                    backend_convert(x, x_tril);

                    auto QP_solver = std::make_shared<utopia::MPGRP<Matrix, Vector> >();

                    Vector lb_tril = local_values(local_size(x_tril).get(0), -9e9);
                    Vector ub_tril = local_values(local_size(x_tril).get(0), 9e9);

                    empty_rhs_tril = 0.0*x_tril;
                    backend_convert(x_eq, x_eq_tril);
                    backend_convert(x_bc_marker, x_bc_marker_tril);

                    QuadraticExtendedFunction<Matrix, Vector> fun_QP_tril(H_tril, g_tril, x_eq_tril, x_bc_marker_tril, empty_rhs_tril);

                    TrustRegionVariableBound<Matrix, Vector> tr_solver(QP_solver);
                    tr_solver.set_box_constraints(make_box_constaints(make_ref(lb_tril), make_ref(ub_tril)));
                    tr_solver.read(input_params_);
                    tr_solver.verbose(verbose_);
                    tr_solver.solve(fun_QP_tril, x_tril);

                #endif //WITH_TRILINOS

            #endif //WITH_PETSC

        }


        template<class Matrix1, class Vector1, class Matrix2, class Vector2>
        std::shared_ptr<QuadraticExtendedFunction<Matrix2, Vector2> > copy_QPfun_to_tril(const Poisson3D<Matrix1, Vector1> & fun)
        {
            Matrix1 H;
            Vector1 g, x_eq, x_bc_marker;

            fun.get_A_rhs(H, g);
            fun.get_eq_constrains_values(x_eq);
            fun.get_eq_constrains_flg(x_bc_marker);

            Matrix2 H_tril;
            Vector2 g_tril, x_eq_tril, x_bc_marker_tril, empty_rhs_tril;
            backend_convert_sparse(H, H_tril);
            backend_convert(g, g_tril);
            backend_convert(x_eq, x_eq_tril);
            backend_convert(x_bc_marker, x_bc_marker_tril);
            empty_rhs_tril = 0.0*x_eq_tril;

            return std::make_shared<QuadraticExtendedFunction<Matrix2, Vector2> >(H_tril, g_tril, x_eq_tril, x_bc_marker_tril, empty_rhs_tril);
        }

        template<class Matrix1, class Vector1, class Matrix2, class Vector2>
        std::shared_ptr<Bratu3D<Matrix2, Vector2> > copy_Bratufun_to_tril(const Poisson3D<Matrix1, Vector1> & fun, const Scalar & lambda=2.1)
        {
            Matrix1 H;
            Vector1 g, x_eq, x_bc_marker;

            fun.get_A_rhs(H, g);
            fun.get_eq_constrains_values(x_eq);
            fun.get_eq_constrains_flg(x_bc_marker);

            Matrix2 H_tril;
            Vector2 x_eq_tril, x_bc_marker_tril, empty_rhs_tril;
            backend_convert_sparse(H, H_tril);
            backend_convert(x_eq, x_eq_tril);
            backend_convert(x_bc_marker, x_bc_marker_tril);
            empty_rhs_tril = 0.0*x_eq_tril;

            int n = pow(fun.dim(), 1./3.);
            double dx = 1./n;

            return std::make_shared<Bratu3D<Matrix2, Vector2> >(H_tril, x_eq_tril, x_bc_marker_tril, empty_rhs_tril, lambda, dx, dx, dx);
        }


        template<class Matrix1, class Vector1, class ProblemType>
        void get_ML_problem(std::vector<std::shared_ptr<Transfer<Matrix1, Vector1> > > & transfers_tril, std::vector<std::shared_ptr<ExtendedFunction<Matrix1, Vector1> > > &  level_functions_tril, Vector1 & x_fine)
        {
            #ifdef WITH_PETSC
                PetscMultilevelTestProblem<PetscMatrix, PetscVector, Poisson3D<PetscMatrix, PetscVector> > multilevel_problem(3, n_levels_, n_);

                auto fun = multilevel_problem.level_functions_[n_levels_-1];
                Poisson3D<PetscMatrix, PetscVector> * fun_Poisson3D = dynamic_cast<Poisson3D<PetscMatrix, PetscVector> *>(fun.get());
                PetscVector x = fun_Poisson3D->initial_guess();

                if(verbose_)
                    fun_Poisson3D->describe();

                backend_convert(x, x_fine);

                for(auto i=0; i < multilevel_problem.transfers_.size(); i++)
                {
                    MatrixTransfer<PetscMatrix, PetscVector> * mat_transfer = dynamic_cast<MatrixTransfer<PetscMatrix, PetscVector> *>(multilevel_problem.transfers_[i].get());

                    Matrix I_tril;
                    backend_convert_sparse(mat_transfer->I(), I_tril);
                    transfers_tril.push_back( std::make_shared<MatrixTransfer<Matrix, Vector> >( std::make_shared<Matrix>(I_tril)));
                }

                level_functions_tril.resize(multilevel_problem.level_functions_.size());
                for(auto i=0; i < multilevel_problem.level_functions_.size(); i++)
                {
                    Poisson3D<PetscMatrix, PetscVector> * fun_Laplace = dynamic_cast<Poisson3D<PetscMatrix, PetscVector> *>(multilevel_problem.level_functions_[i].get());

                    if(std::is_same<ProblemType, Bratu3D<Matrix1, Vector1>>::value)
                    {
                        level_functions_tril[i] = copy_Bratufun_to_tril<PetscMatrix, PetscVector, Matrix, Vector>(*fun_Laplace);
                    }
                    else
                    {
                        level_functions_tril[i] = copy_QPfun_to_tril<PetscMatrix, PetscVector, Matrix, Vector>(*fun_Laplace);
                    }
                }

            #endif //WITH_PETSC

        }



        void RMTR_inf_test()
        {
            Vector x_fine;
            std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > transfers_tril;
            std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >  level_functions_tril;

            get_ML_problem<Matrix, Vector, Bratu3D<Matrix, Vector>>(transfers_tril, level_functions_tril, x_fine);
            // get_ML_problem<Matrix, Vector, Poisson3D<Matrix, Vector>>(transfers_tril, level_functions_tril, x_fine);

            auto tr_strategy_coarse = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
            auto tr_strategy_fine   = std::make_shared<utopia::MPGRP<Matrix, Vector> >();

            auto rmtr = std::make_shared<RMTR_inf<Matrix, Vector, GALERKIN> >(n_levels_);

            // Set TR-QP strategies
            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);

            // Transfers and objective functions
            rmtr->set_transfer_operators(transfers_tril);
            rmtr->set_functions(level_functions_tril);

            // add constraints
            Vector lb = local_values(local_size(x_fine), -9e9);
            Vector ub = local_values(local_size(x_fine), 1.1);

            auto box = make_box_constaints(make_ref(lb), make_ref(ub));
            rmtr->set_box_constraints(box);


            rmtr->read(input_params_);
            rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);
            rmtr->verbose(verbose_);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

            // Solve
            rmtr->solve(x_fine);
        }


        void RMTR_l2_test()
        {
            Vector x_fine;
            std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > transfers_tril;
            std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >  level_functions_tril;

            // get_ML_problem<Matrix, Vector, Bratu3D<Matrix, Vector>>(transfers_tril, level_functions_tril, x_fine);
            get_ML_problem<Matrix, Vector, Poisson3D<Matrix, Vector>>(transfers_tril, level_functions_tril, x_fine);


            auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
            auto precond = std::make_shared<GaussSeidel<Matrix, Vector, HOMEMADE> >();
            precond->max_it(1);
            tr_strategy_fine->set_preconditioner(precond);
            // tr_strategy_fine->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());


            auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
            // auto precond_coarse = std::make_shared<GaussSeidel<Matrix, Vector, HOMEMADE> >();
            // precond_coarse->max_it(1);
            // tr_strategy_coarse->set_preconditioner(precond);
            tr_strategy_coarse->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());


            auto rmtr = std::make_shared<RMTR<Matrix, Vector, GALERKIN> >(n_levels_);

            // Set TR-QP strategies
            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);

            // Transfers and objective functions
            rmtr->set_transfer_operators(transfers_tril);
            rmtr->set_functions(level_functions_tril);

            rmtr->read(input_params_);
            rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);
            rmtr->verbose(verbose_);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

            // Solve
            rmtr->solve(x_fine);
        }


        void Quasi_RMTR_l2_test()
        {
            Vector x_fine;
            std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > transfers_tril;
            std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >  level_functions_tril;
            get_ML_problem<Matrix, Vector, Poisson3D<Matrix, Vector>>(transfers_tril, level_functions_tril, x_fine);

            auto rmtr = std::make_shared<QuasiRMTR<Matrix, Vector, FIRST_ORDER> >(n_levels_);

            // auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
            // tr_strategy_fine->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());

            // auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
            // tr_strategy_coarse->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());

            // Set TR-QP strategies
            // rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            // rmtr->set_fine_tr_strategy(tr_strategy_fine);

            // // here, we can set also hessian approx strategy for individually
            // SizeType memory_size = 5;
            // auto hes_approx   = std::make_shared< LBFGS<Vector> >(memory_size);
            // rmtr->set_hessian_approximation_strategy();

            const SizeType memory_size = 5;
            std::vector<std::shared_ptr<HessianApproximation<Vector> > > hess_approxs(n_levels_);
            for(auto l=0; l < n_levels_; l++)
            {
                auto hes_approx   = std::make_shared<LBFGS<Vector> >(memory_size);
                hess_approxs[l] = hes_approx;
            }
            rmtr->set_hessian_approximation_strategies(hess_approxs);


            std::vector<std::shared_ptr<utopia::MatrixFreeTRSubproblem<Vector> > > subproblems(n_levels_);
            for(auto l=0; l < n_levels_; l++)
            {
                auto tr_strategy = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
                tr_strategy->atol(1e-14);

                auto precond = hess_approxs[l]->build_Hinv_precond();
                tr_strategy->set_preconditioner(precond);

                subproblems[l] = tr_strategy;
            }
            rmtr->set_tr_strategies(subproblems);

            // Transfers and objective functions
            rmtr->set_transfer_operators(transfers_tril);
            rmtr->set_functions(level_functions_tril);

            rmtr->read(input_params_);
            rmtr->pre_smoothing_steps(10);
            rmtr->post_smoothing_steps(10);
            rmtr->max_coarse_it(10);

            rmtr->max_sucessful_smoothing_it(5);
            rmtr->max_sucessful_coarse_it(10);


            rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);
            rmtr->verbose(verbose_);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

            // Solve
            rmtr->solve(x_fine);
        }


        void Quasi_RMTR_inf_test()
        {
            Vector x_fine;
            std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > transfers_tril;
            std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >  level_functions_tril;
            get_ML_problem<Matrix, Vector, Poisson3D<Matrix, Vector>>(transfers_tril, level_functions_tril, x_fine);

            auto rmtr = std::make_shared<QuasiRMTR_inf<Matrix, Vector, FIRST_ORDER> >(n_levels_);

            auto tr_strategy_fine = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
            tr_strategy_fine->atol(1e-12);
           //  tr_strategy_fine->verbose(true);
            auto tr_strategy_coarse = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
            tr_strategy_coarse->atol(1e-12);

            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);

            // here, we can set also hessian approx strategy for individually
            SizeType memory_size = 5;
            auto hes_approx   = std::make_shared< LBFGS<Vector> >(memory_size);
            rmtr->set_hessian_approximation_strategy(hes_approx);

            // Transfers and objective functions
            rmtr->set_transfer_operators(transfers_tril);
            rmtr->set_functions(level_functions_tril);

            // add constraints
            // Vector lb = local_values(local_size(x_fine), -9e9);
            // Vector ub = local_values(local_size(x_fine), 1.1);

            // auto box = make_box_constaints(make_ref(lb), make_ref(ub));
            // rmtr->set_box_constraints(box);

            rmtr->read(input_params_);
            rmtr->pre_smoothing_steps(10);
            rmtr->post_smoothing_steps(10);
            rmtr->max_coarse_it(10);

            rmtr->max_sucessful_smoothing_it(5);
            rmtr->max_sucessful_coarse_it(10);

            rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);
            rmtr->verbose(verbose_);
           //  rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

            // Solve
            rmtr->solve(x_fine);

            // disp(x_fine);
        }

        void for_each_loop_test()
        {
            Vector x = values(n_, 2.);
            Vector y = values(n_, 1.);
            Vector z = values(n_, 0.);

            using ForLoop = utopia::ParallelFor<Traits<Vector>::Backend>;

           {
                auto d_x = const_device_view(x);
                auto d_y = const_device_view(y);
                auto d_z = device_view(z);

                ForLoop::apply(range(z), UTOPIA_LAMBDA(const SizeType i)
                {
                    const Scalar xi = d_x.get(i);
                    const Scalar yi = d_y.get(i);
                    d_z.set(i, xi - yi);
                });
            }

            Scalar sum_z = sum(z);
            utopia_test_assert(approxeq(Scalar(n_), sum_z));
        }

        void parallel_each_write_test()
        {
            Vector x = values(n_, 2.);
            Vector y = values(n_, 1.);
            Vector z = values(n_, 0.);

            {
                auto d_x = const_device_view(x);
                auto d_y = const_device_view(y);
                auto d_z = device_view(z);

                parallel_each_write(z, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    const Scalar xi = d_x.get(i);
                    const Scalar yi = d_y.get(i);
                    return xi - yi;
                });
            }

            Scalar sum_z = sum(z);
            utopia_test_assert(approxeq(Scalar(n_), sum_z));
        }
        
        void multi_reduce_test()
        {
            Vector x = values(n_, 2.);
            Vector y = values(n_, 1.);

            each_write(x, [](const SizeType &i) -> Scalar {
                return -(i + 1.0);
            });

            each_write(y, [](const SizeType &i) -> Scalar {
                return (i + 1.0);
            });

            const Scalar m = multi_min(x, y);
            utopia_test_assert(approxeq(m, Scalar(-n_)));
        }

    private:
        SizeType n_;
        SizeType n_levels_;
        Scalar lambda_;
        bool verbose_;
        bool output_vtk_;

        InputParameters input_params_;

    };

    void hck()
    {
#ifdef WITH_PETSC
        auto n_levels    = 3;

        auto coarse_dofs = 1000;
        auto verbose     = true;

        // HckTests<PetscMatrix, PetscVector>(coarse_dofs, n_levels, 1.0, false, true).run_petsc();
        // HckTests<PetscMatrix, PetscVector>(coarse_dofs, n_levels, 1.0, verbose, true).run_trilinos();

#ifdef WITH_TRILINOS
        HckTests<TpetraMatrixd, TpetraVectord>(coarse_dofs, n_levels, 1.0, verbose, true).run_trilinos();
#endif
#endif

    }

    UTOPIA_REGISTER_TEST_FUNCTION(hck);
}

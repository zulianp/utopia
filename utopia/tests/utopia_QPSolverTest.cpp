#include "utopia_Testing.hpp"

#include "utopia.hpp"

#include "utopia_ProjectedConjugateGradient.hpp"
#include "utopia_ProjectedGradient.hpp"

#include "test_problems/utopia_QPSolverTestProblem.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

#include "utopia_polymorphic_QPSolver.hpp"

#include "utopia_MonotoneAlgebraicMultigrid.hpp"
#include "utopia_MonotoneMultigrid.hpp"
#include "utopia_MultigridQR.hpp"

#include "utopia_MultilevelTestProblem1D.hpp"
#include "utopia_Poisson1D.hpp"

#include "utopia_Agglomerate.hpp"
#include "utopia_LogBarrierFunction.hpp"
#include "utopia_LogBarrierQPSolver.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"
#endif  // UTOPIA_WITH_PETSC

namespace utopia {

    template <class Matrix, class Vector>
    class QPSolverTest {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;
        using IndexArray = typename Traits::IndexArray;
        using IndexSet = typename Traits::IndexSet;

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

        template <class QPSolver>
        void run_qp_solver(QPSolver &qp_solver) const {
            QPSolverTestProblem<Matrix, Vector>::run(n, verbose, qp_solver);
        }

        void pg_test() const {
            ProjectedGradient<Matrix, Vector> pg;
            run_qp_solver(pg);
        }

        void pcg_test() const {
            ProjectedConjugateGradient<Matrix, Vector> pcg;
            run_qp_solver(pcg);
        }

        void ngs_test() const {
            ProjectedGaussSeidel<Matrix, Vector> pgs;
            run_qp_solver(pgs);
        }

        void nblockgs_test() const {
            InputParameters params;
            params.set("block_size", 2);

            ProjectedGaussSeidel<Matrix, Vector> pgs;
            pgs.read(params);

            run_qp_solver(pgs);
        }

        static void create_symm_lapl_test_data(Comm &comm,
                                               Matrix &A,
                                               Vector &b,
                                               BoxConstraints<Vector> &box,
                                               SizeType n = 100,
                                               const bool boundary_conds = true) {
            A.sparse(layout(comm, Traits::decide(), Traits::decide(), n, n), 3, 2);
            assemble_symmetric_laplacian_1D(A, boundary_conds);

            auto h = 1. / (n - 1.);
            A = 1. / h * A;

            {
                Range r = row_range(A);
                const SizeType r_begin = r.begin();
                const SizeType r_end = r.end();

                Write<Matrix> w(A);
                if (r_begin == SizeType(0)) {
                    A.set(0, 0, 1.);
                }

                if (r_end == n) {
                    A.set(n - 1, n - 1, 1.);
                }
            }

            b.values(row_layout(A), 50.0);

            {
                Range row_range = range(b);
                const SizeType r_begin = row_range.begin();
                const SizeType r_end = row_range.end();

                Write<Vector> w(b);

                for (SizeType r = r_begin; r != r_end; ++r) {
                    if (r >= n / 2.) {
                        b.set(r, -50.0);
                    }
                    if (r == 0) {
                        b.set(r, 0);
                    }

                    if (r == (n - 1)) {
                        b.set(r, 0);
                    }
                }
            }

            b = h * b;

            auto lb = std::make_shared<Vector>(row_layout(A), -0.5);
            auto ub = std::make_shared<Vector>(row_layout(A), 0.5);

            box = make_box_constaints(lb, ub);
        }

        void log_barrier_test() {
            auto &&comm = Comm::get_default();

            Matrix A;
            Vector b;
            BoxConstraints<Vector> box;
            create_symm_lapl_test_data(comm, A, b, box);

            // b *= 0.5;
            // box.lower_bound() = nullptr;
            // box.upper_bound() = nullptr;
            QuadraticFunction<Matrix, Vector> fun(make_ref(A), make_ref(b));

            InputParameters params;
            // params.set("verbose", true);
            params.set("barrier_parameter", 1e-5);
            params.set("barrier_parameter_shrinking_factor", 0.7);

            LogBarrierFunction<Matrix, Vector> barrier(make_ref(fun), make_ref(box));
            barrier.read(params);

            ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector>>());
            // cg.verbose(true);
            cg.max_it(20);

            Newton<Matrix, Vector> newton(make_ref(cg));
            // newton.verbose(true);

            Vector x(layout(b));

            // Linear solve first to get closer to solution
            cg.solve(A, b, x);
            barrier.project_onto_feasibile_region(x);
            newton.solve(barrier, x);
        }

        void log_barrier_qp_solver_test(const std::string &barrier_function_type, const bool verbose) {
            auto &&comm = Comm::get_default();

            Matrix A;
            Vector b;
            BoxConstraints<Vector> box;
            create_symm_lapl_test_data(comm, A, b, box);
            box.lower_bound() = nullptr;

            InputParameters params;
            // params.set("verbose", true);
            params.set("barrier_parameter", 1e-5);
            params.set("barrier_thickness", 0.01);
            params.set("barrier_parameter_shrinking_factor", 0.7);
            params.set("min_barrier_parameter", 1e-8);
            params.set("verbose", verbose);
            params.set("function_type", barrier_function_type);
            params.set("max-it", 10);

            LogBarrierQPSolver<Matrix, Vector> solver;
            solver.set_box_constraints(box);
            solver.read(params);

            Vector x(layout(b), 0.);
            utopia_test_assert(solver.solve(A, b, x));

            if (Traits::Backend == PETSC) {
                rename("x", x);
                write("X.m", x);

                if (box.has_lower_bound()) {
                    rename("lb", *box.lower_bound());
                    write("LB.m", *box.lower_bound());
                }

                if (box.has_upper_bound()) {
                    rename("ub", *box.upper_bound());
                    write("UB.m", *box.upper_bound());
                }

                rename("a", A);
                write("A.m", A);

                rename("b", b);
                write("B.m", b);
            }
        }

        void log_barrier_qp_solver_test() {
            // log_barrier_qp_solver_test("LogBarrierFunction", true);
            log_barrier_qp_solver_test("LogBarrierFunction", false);
            // log_barrier_qp_solver_test("BoundedLogBarrierFunction", true);

            // Utopia::Abort("BYE");
        }

        void interior_point_qp_solver_test() {
            auto &&comm = Comm::get_default();

            SizeType n = 2000;
            if (Traits::Backend == BLAS) {
                n = 20;
            }

            Matrix A, B;
            Vector b;
            BoxConstraints<Vector> box;
            create_symm_lapl_test_data(comm, A, b, box, n);

            Vector x(layout(b), 0.), lambda(layout(b), 0.), slack(layout(b), 0.);
            Vector buff(layout(b), 0.);
            Scalar sigma = 0;  // in [0, 1]
            Scalar mu = 0.0;

            B.identity(layout(A), 1.0);

            ConjugateGradient<Matrix, Vector, HOMEMADE> solver;
            solver.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector>>());
            // solver.verbose(true);

            auto duality_measure = [&](const Vector &lambda, const Vector &s) -> Scalar {
                return dot(lambda, s) / size(lambda);
            };

            auto residual_d = [&](const Vector &x, const Vector &lambda, Vector &r_d) -> bool {
                r_d = A * x + transpose(B) * lambda - b;
                return true;
            };

            auto residual_p = [&](const Vector &x, const Vector &s, Vector &r_p) -> bool {
                r_p = B * x - *box.upper_bound() + s;
                return true;
            };

            auto residual_s = [&](const Vector &lambda, const Vector &s, Vector &r_s) -> bool {
                buff.set(sigma * mu);
                r_s = e_mul(s, lambda);
                r_s -= buff;
                return true;
            };

            auto grad_f = [&A, &b](const Vector &x, Vector &g) -> bool {
                g = A * x - b;
                return true;
            };

            auto hessian_f = [&A](const Vector &x, Matrix &H) -> bool {
                H = A;
                return true;
            };

            auto is_positive = [](const Vector &x) -> bool {
                auto x_view = local_view_device(x);

                int n_violations = 0;
                parallel_reduce(
                    local_range_device(x),
                    UTOPIA_LAMBDA(const SizeType i)->int {
                        auto xi = x_view.get(i);
                        // if (xi <= 0) {
                        //     std::cout << "x(" << i << ") == " << xi << "\n";
                        // }

                        return (xi <= 0) ? 1 : 0;
                    },
                    n_violations);

                n_violations = x.comm().sum(n_violations);

                // if (n_violations) utopia::out() << "n_violations: " << n_violations << "/" << x.size() << "\n";

                return n_violations == 0;
            };

            auto apply_bc = [](Vector &x) {
                auto r = range(x);
                auto N = x.size();

                Write<Vector> w(x);
                if (r.inside(0)) {
                    x.set(0, 0);
                }

                if (r.inside(N - 1)) {
                    x.set(N - 1, 0);
                }
            };

            auto compute_alpha_with_tau = [&buff](const Vector &x, const Vector &delta_x, const Scalar tau) -> Scalar {
                {
                    auto x_view = local_view_device(x);
                    auto delta_x_view = local_view_device(delta_x);

                    buff.set(1.);

                    auto buff_view = local_view_device(buff);

                    parallel_for(
                        local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                            auto xi = x_view.get(i);
                            auto dxi = delta_x_view.get(i);

                            if (xi + dxi < (1 - tau) * xi) {
                                auto val = -xi * tau / dxi;
                                buff_view.set(i, val);
                            }
                        });
                }

                return min(buff);
            };

            auto compute_alpha = [&buff](const Vector &x, const Vector &delta_x) -> Scalar {
                {
                    auto x_view = local_view_device(x);
                    auto delta_x_view = local_view_device(delta_x);

                    buff.set(1.);

                    auto buff_view = local_view_device(buff);

                    parallel_for(
                        local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                            auto xi = x_view.get(i);
                            auto dxi = delta_x_view.get(i);

                            if (dxi < 0 && xi != 0) {
                                auto val = -xi / dxi;
                                buff_view.set(i, val);
                            }
                        });
                }

                return min(buff);
            };

            Matrix H;
            Vector g(layout(b), 0.), delta_x(layout(b), 0.), delta_lambda, delta_slack, rd, rp, rs, rs_aff, slack_inv,
                lambda_inv;

            x = min(x, *box.upper_bound());
            lambda.values(layout(x), 1.);

            slack = *box.upper_bound() - B * x;
            slack.e_max(1e-8);  // 0.1

            assert(is_positive(slack));
            assert(is_positive(lambda));

            bool verbose = Traits::Backend == PETSC;
            int max_it = 1000;
            Scalar tol = 1e-8;
            for (int k = 0; k < max_it; ++k) {
                mu = duality_measure(lambda, slack);
                sigma = 0;

                residual_d(x, lambda, rd);
                residual_p(x, slack, rp);
                residual_s(lambda, slack, rs);

                H = A;
                e_pseudo_inv(slack, slack_inv);
                e_pseudo_inv(lambda, lambda_inv);

                /// FIXME begin
                Vector temp_vec = e_mul(slack_inv, lambda);
                apply_bc(temp_vec);

                Matrix temp_mat = diag(temp_vec);

                Matrix temp_mat_2 = transpose(B) * temp_mat * B;
                H += temp_mat_2;
                /// FIXME end

                g = -transpose(B) * e_mul(e_mul(rp, lambda) - rs, slack_inv);
                g -= rd;

                apply_bc(g);

                delta_x.set(0.);
                if (!solver.solve(H, g, delta_x)) {
                    // break;
                }

                //////////////////////////////////////////////////////////////////////

                delta_lambda = e_mul(slack_inv, e_mul(lambda, B * delta_x) + e_mul(lambda, rp) - rs);
                delta_slack = -e_mul(lambda_inv, rs + e_mul(slack, delta_lambda));

                // /////////////////////////////////////////////////////////////////////////
                // // https://en.wikipedia.org/wiki/Mehrotra_predictor%E2%80%93corrector_method

                Scalar alpha_primal = compute_alpha(slack, delta_slack);
                assert(alpha_primal == alpha_primal);

                Scalar alpha_dual = compute_alpha(lambda, delta_lambda);
                assert(alpha_dual == alpha_dual);

                Scalar alpha = std::min(alpha_primal, alpha_dual);

                // Scalar mu_aff = dot(slack + alpha_primal * delta_slack, lambda + alpha_dual * delta_lambda);
                Scalar mu_aff = dot(slack + alpha * delta_slack, lambda + alpha * delta_lambda);
                mu_aff /= size(lambda);
                sigma = pow(mu_aff / mu, 3);
                assert(sigma == sigma);

                // //////////////////////////////////////

                buff.set(sigma * mu);
                rs_aff = e_mul(slack, lambda) + e_mul(delta_slack, delta_lambda) - buff;

                g = -transpose(B) * e_mul(slack_inv, e_mul(lambda, rp) - rs_aff);
                g -= rd;

                apply_bc(g);

                delta_x.set(0.);
                if (!solver.solve(H, g, delta_x)) {
                    // break;
                }

                delta_lambda = e_mul(slack_inv, e_mul(lambda, B * delta_x) + e_mul(lambda, rp) - rs_aff);
                delta_slack = -e_mul(lambda_inv, rs_aff + e_mul(slack, delta_lambda));

                // ////////////////////////////////////////
                Scalar tau = 1;
                alpha_primal = compute_alpha_with_tau(slack, delta_slack, tau);
                assert(alpha_primal == alpha_primal);

                alpha_dual = compute_alpha_with_tau(lambda, delta_lambda, tau);
                assert(alpha_dual == alpha_dual);
                alpha = std::min(alpha_primal, alpha_dual);

                x += alpha * delta_x;
                lambda += alpha * delta_lambda;
                slack += alpha * delta_slack;

                Scalar norm_dx = norm1(delta_x);

                if (verbose) {
                    std::stringstream ss;
                    ss << "==============================\n";
                    ss << "iteration:\t" << k << "\n";
                    ss << "norm_dx:\t" << norm_dx << "\n";
                    ss << "alpha:\t" << alpha << "\n";
                    ss << "alpha_primal:\t" << alpha_primal << "\n";
                    ss << "alpha_dual:\t" << alpha_dual << "\n";
                    ss << "mu:\t" << mu << "\n";
                    ss << "sigma:\t" << sigma << "\n";

                    x.comm().root_print(ss.str());
                }

                bool slack_is_positive = is_positive(slack);

                if (((k == (max_it - 1)) || norm_dx < tol)) {
                    if (norm_dx < tol && verbose) {
                        std::stringstream ss;
                        ss << "Solution found in " << k << " iterations\n";
                        x.comm().root_print(ss.str());
                    }

                    if (Traits::Backend == PETSC && verbose) {
                        rename("x", x);
                        rename("lambda", lambda);
                        rename("slack", slack);

                        write("IP_X.m", x);
                        write("IP_S.m", slack);
                        write("IP_L.m", lambda);
                    }

                    break;
                }
            }
        }

        void MPRGP_test() const {
            MPRGP<Matrix, Vector> qp_solver;
            run_qp_solver(qp_solver);

            auto &&comm = Comm::get_default();

            Matrix A;
            Vector b;
            BoxConstraints<Vector> box;
            create_symm_lapl_test_data(comm, A, b, box);

            Vector x = 0 * b;

            qp_solver.set_box_constraints(box);
            qp_solver.verbose(false);
            qp_solver.max_it(n * 2);
            qp_solver.set_eig_comp_tol(1e-1);
            qp_solver.solve(A, b, x);
        }

        void MG_QR_test() {
            bool verbose = false;

            Vector rhs, x;
            Vector upper_bound, lower_bound;
            Matrix A, R, Q, Ih_fine, Rot;
            Matrix Ih1, Ih0;

            const std::string data_path = Utopia::instance().get("data_path");

            read(data_path + "/forQR/b", rhs);
            read(data_path + "/forQR/x", x);
            read(data_path + "/forQR/A", A);
            read(data_path + "/forQR/Q", Q);
            read(data_path + "/forQR/R", R);
            read(data_path + "/forQR/Rot", Rot);
            read(data_path + "/forQR/ub", upper_bound);
            read(data_path + "/forQR/lb", lower_bound);

            read(data_path + "/forQR/Ih", Ih_fine);
            read(data_path + "/forQR/I2h", Ih1);
            read(data_path + "/forQR/I3h", Ih0);

            x.set(1);

            auto num_levels = 3;

            // chop_abs(Q, 1e-7);

            // std::cout<<"A: "<< local_size(A).get(0) << "  \n";
            // std::cout<<"Q: "<< local_size(Q).get(0) << "  \n";
            // std::cout<<"R: "<< local_size(R).get(0) << ","<<local_size(R).get(1) <<
            // "  \n"; std::cout<<"I: "<< local_size(Ih_fine).get(0) << "  \n";
            // std::cout<<"rhs: "<< local_size(rhs).get(0) << "  \n";

            R = transpose(R);

            // version 1
            Matrix QtAQ = transpose(Q) * Rot * A * Rot * Q;
            Matrix QtIh = transpose(Q) * Rot * Ih_fine;
            Vector Qtrhs = transpose(Q) * Rot * rhs;
            Vector Qtx = transpose(Q) * Rot * x;

            // Matrix QtAQ  = Rot*A*Rot;
            // Matrix QtIh  = Rot* Ih_fine;
            // Vector Qtrhs = Rot *rhs;
            // Vector Qtx   = Rot *x;

            auto fine_smoother = std::make_shared<ProjectedGaussSeidelQR<Matrix, Vector>>();
            fine_smoother->set_R(R);  // Monotone

            auto coarse_smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();
            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>("mumps", "lu");
            // MultigridQR<Matrix, Vector> multigrid(fine_smoother, coarse_smoother, direct_solver, num_levels);  // QR
            MonotoneMultigrid<Matrix, Vector> multigrid(fine_smoother, coarse_smoother, direct_solver);  // Monotone

            std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> interpolation_operators;
            interpolation_operators.resize(num_levels - 1);
            interpolation_operators[1] =
                std::make_shared<IPTruncatedTransfer<Matrix, Vector>>(std::make_shared<Matrix>(QtIh));
            interpolation_operators[0] = std::make_shared<IPTransfer<Matrix, Vector>>(std::make_shared<Matrix>(Ih1));

            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.max_it(40);
            multigrid.pre_smoothing_steps(3);
            multigrid.post_smoothing_steps(3);
            multigrid.verbose(verbose);
            // multigrid.verbose(true);

            // multigrid.mg_type(2);

            // This should be somewhere else...
            // multigrid.set_QR(Q, R);  // QR
            multigrid.set_box_constraints(make_box_constaints(make_ref(lower_bound), make_ref(upper_bound)));

            multigrid.solve(QtAQ, Qtrhs, Qtx);
            x = Rot * Q * Qtx;
            // x = Rot*Qtx;

            // write("x.m", x);
            // write("IX.m", Qtx);
            // disp(x);
        }

        void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(interior_point_qp_solver_test);
            UTOPIA_RUN_TEST(log_barrier_test);
            UTOPIA_RUN_TEST(log_barrier_qp_solver_test);
            UTOPIA_RUN_TEST(pg_test);
            UTOPIA_RUN_TEST(pcg_test);
            UTOPIA_RUN_TEST(ngs_test);
            UTOPIA_RUN_TEST(MPRGP_test);
            UTOPIA_RUN_TEST(nblockgs_test);
        }

        void run_GS_QR() {
            if (mpi_world_size() > 1) return;

            print_backend_info();
            UTOPIA_RUN_TEST(MG_QR_test);
        }

        QPSolverTest() : n(20) {}

        SizeType n = 20;
        bool verbose = false;
    };

    template <class Matrix, class Vector>
    class BDDOperator : public Operator<Vector> {
    public:
        using Communicator = typename Traits<Vector>::Communicator;

        static void parallel_to_serial(const Vector &x_from, Vector &x_to) {
            {
                // FIXME (copying stuff just because of abstractions!)
                auto x_from_view = local_view_device(x_from);
                auto x_to_view = local_view_device(x_to);
                parallel_for(
                    local_range_device(x_to),
                    UTOPIA_LAMBDA(const SizeType i) { x_to_view.set(i, x_from_view.get(i)); });
            }
        }

        static void serial_to_parallel(const Vector &x_from, Vector &x_to) { parallel_to_serial(x_from, x_to); }

        bool apply(const Vector &x_G, Vector &rhs_G) const override {
            parallel_to_serial(x_G, *xL_);

            *A_IG_x_ = (*A_IG_) * (*xL_);
            sol_I_->set(0);

            A_II_inv_->apply(*A_IG_x_, *sol_I_);

            (*rhsL_) = (*A_GI_) * (*sol_I_);

            if (empty(rhs_G)) {
                rhs_G.zeros(layout(x_G));
            }

            (*rhsL_) *= -1;

            serial_to_parallel(*rhsL_, rhs_G);
            rhs_G += (*A_GG_) * x_G;
            return false;
        }

        Size size() const override { return A_GG_->size(); }

        Size local_size() const override { return A_GG_->local_size(); }

        Communicator &comm() override { return A_GG_->comm(); }
        const Communicator &comm() const override { return A_GG_->comm(); }

        void init(const std::shared_ptr<Matrix> &A_GG,
                  const std::shared_ptr<Matrix> &A_GI,
                  const std::shared_ptr<Matrix> &A_II,
                  const std::shared_ptr<Matrix> &A_IG) {
            A_GG_ = A_GG;
            A_GI_ = A_GI;
            A_II_ = A_II;
            A_IG_ = A_IG;

            A_II_inv_ = std::make_shared<Factorization<Matrix, Vector>>("superlu", "lu");
            A_II_inv_->update(A_II_);
        }

        void init_rhs(const Vector &rhs_G, const Vector &rhs_I) {
            assert(A_II_inv_);

            // Compute rhs
            secant_G_ = std::make_shared<Vector>(layout(rhs_G));
            xL_ = std::make_shared<Vector>(row_layout(*A_GI_));
            rhsL_ = std::make_shared<Vector>(row_layout(*A_GI_));
            A_IG_x_ = std::make_shared<Vector>(row_layout(*A_II_));
            sol_I_ = std::make_shared<Vector>(row_layout(*A_II_));

            Vector inv_A_II_rhs_I(layout(rhs_I), 0);
            A_II_inv_->apply(rhs_I, inv_A_II_rhs_I);

            (*secant_G_) = rhs_G;

            Vector temp = (*A_GI_) * inv_A_II_rhs_I;

            assert(secant_G_->local_size() == temp.local_size());

            {
                auto sG_view = local_view_device(*secant_G_);
                auto temp_view = local_view_device(temp);

                parallel_for(
                    local_range_device(temp),
                    UTOPIA_LAMBDA(const SizeType i) { sG_view.set(i, sG_view.get(i) - temp_view.get(i)); });
            }

            // disp(rhs_G);
            // disp(*secant_G_);

            // for (SizeType r = 0; r < comm().size(); ++r) {
            //     comm().barrier();

            //     if (r == comm().rank()) {
            //         std::cout << "=========================================\n";
            //         std::cout << "Rank: " << r << "\n";
            //         std::cout << "=========================================\n";
            //         std::cout << "A_GI * A_II^-1 * rhs_I:\n";
            //         disp(temp);

            //         std::cout << "\n" << std::flush;
            //     }
            // }
        }

        void finalize(const Vector &x_G, const Vector &rhs_I, Vector &x_I) {
            parallel_to_serial(x_G, *xL_);

            *A_IG_x_ = (*A_IG_) * (*xL_);
            *rhsL_ = rhs_I;
            *rhsL_ -= *A_IG_x_;

            if (empty(x_I)) {
                x_I.zeros(layout(*rhsL_));
            } else {
                x_I.set(0);
            }

            A_II_inv_->apply(*rhsL_, x_I);
        }

        std::shared_ptr<Matrix> A_GG_, A_GI_, A_II_, A_IG_;
        std::shared_ptr<Vector> secant_G_;
        std::shared_ptr<Factorization<Matrix, Vector>> A_II_inv_;
        std::shared_ptr<Vector> xL_, rhsL_, A_IG_x_, sol_I_;
    };

    // FIXME merge with the other once it is poperly implemented
    template <class Matrix, class Vector>
    class PQPSolverTest {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;
        using IndexArray = typename Traits::IndexArray;
        using IndexSet = typename Traits::IndexSet;

        void run() {
            UTOPIA_RUN_TEST(MPRGP_DD);
            // FIXME
            UTOPIA_RUN_TEST(poly_qp);
        }

        void MPRGP_DD() {
            auto &&comm = Comm::get_default();

            static const bool verbose = false;

            Matrix A;
            Vector b;
            BoxConstraints<Vector> box;

            Vector oracle;

            std::stringstream c_ss;
            Chrono c;

            if (true) {
                c.start();

                SizeType n = 1e3;
                QPSolverTest<Matrix, Vector>::create_symm_lapl_test_data(comm, A, b, box, n, true);

                c.stop();
                c_ss << "Problem initialization\n" << c << "\n";

                if (n <= 1e5) {
                    c.start();
                    Factorization<Matrix, Vector> solver;
                    solver.solve(A, b, oracle);
                    c.stop();
                    c_ss << "Direct solver\n" << c << "\n";
                }

            } else {
                c.start();

                Path dir = "../data/test/CG_DD/mats_tests_2d_tri3";
                read(dir / "A", A);
                read(dir / "b", b);
                read(dir / "x", oracle);

                c.stop();
                c_ss << "Read problem from disk\n" << c << "\n";
            }

            c.start();

            // disp(A);

            Vector x(layout(b), 0.);

            auto r = row_range(A);
            SizeType n_local = r.extent();
            IndexSet min_idx(n_local, A.rows()), max_idx(n_local, 0);

            A.read([&](const SizeType i, const SizeType j, const Scalar) {
                min_idx[i - r.begin()] = std::min(SizeType(min_idx[i - r.begin()]), j);
                max_idx[i - r.begin()] = std::max(SizeType(max_idx[i - r.begin()]), j);
            });

            SizeType n_interface = 0;
            // std::stringstream ss;
            for (SizeType i = 0; i < n_local; ++i) {
                if (!r.inside(min_idx[i]) || !r.inside(max_idx[i])) {
                    ++n_interface;
                    // ss << (i + r.begin()) << " ";
                }
            }

            IndexSet local_interface_idx(n_interface, 0);
            IndexSet global_to_local(n_local, -1);
            IndexSet local_to_global(n_interface, 0);

            n_interface = 0;
            for (SizeType i = 0; i < n_local; ++i) {
                if (!r.inside(min_idx[i]) || !r.inside(max_idx[i])) {
                    local_to_global[n_interface] = (i + r.begin());
                    global_to_local[i] = n_interface;
                    local_interface_idx[n_interface++] = i;
                }
            }

            SizeType interface_offset = 0;
            x.comm().exscan_sum(&n_interface, &interface_offset, 1);

            // x.comm().synched_print(ss.str());

            Matrix A_II_view;
            local_block_view(A, A_II_view);

            Matrix A_II = A_II_view;
            set_zero_rows(A_II, local_interface_idx, 1);

            A_II.transform_ijv([&](const SizeType i, const SizeType j, const Scalar value) -> Scalar {
                bool is_interface = global_to_local[j] != -1;
                if (is_interface) {
                    return i == j ? 1.0 : 0.;
                } else {
                    return value;
                }
            });

            // for (SizeType r = 0; r < comm.size(); ++r) {
            //     comm.barrier();

            //     if (r == comm.rank()) {
            //         std::cout << "interface_offset: " << interface_offset << std::endl;

            //         for (SizeType i = 0; i < n_interface; ++i) {
            //             std::cout << (interface_offset + i) << " -> " << local_to_global[i] << " ";
            //         }

            //         std::cout << "\n";
            //     }
            // }

            // comm.barrier();

            const PetscInt *ranges;
            MatGetOwnershipRanges(A.raw_type(), &ranges);

            SizeType n_col_dofs = A.cols();
            int comm_size = comm.size();
            auto find_rank = [=](const SizeType i) -> int {
                int rank = std::min(int(i * (float(comm_size) / n_col_dofs)), comm_size - 1);

                bool found = (i >= ranges[rank]) && (i < ranges[rank + 1]);

                while (!found) {
                    if (i < ranges[rank]) {
                        --rank;
                    } else if (i >= ranges[rank + 1]) {
                        ++rank;
                    } else {
                        assert(false);
                    }

                    // if (i >= ranges[comm_size]) {
                    //     std::cout << comm.rank() << ": [" << ranges[0] << "-" << ranges[comm.size()] << ")"
                    //               << " r: " << rank << " " << i << "\n";
                    // }

                    assert(i >= 0);
                    assert(i < ranges[comm_size]);
                    assert(rank < comm_size);
                    assert(rank >= 0);

                    found = (i >= ranges[rank]) && (i < ranges[rank + 1]);
                }

                assert(found);

                return rank;
            };

            std::vector<SizeType> counter(comm.size(), 0);
            std::vector<IndexSet> recv_list(comm.size()), send_list(comm.size());
            std::vector<SizeType> interface_ghosts;

            {
                Mat d, o;
                MatMPIAIJGetSeqAIJ(A.raw_type(), &d, &o, nullptr);

                utopia::PetscCrsView d_crs_view(d);
                utopia::PetscCrsView o_crs_view(o);

                PetscInt nghosts;
                const PetscInt *ghosts = nullptr;
                MatGetGhosts(A.raw_type(), &nghosts, &ghosts);
                interface_ghosts.resize(nghosts, -1);

                auto d_row_ptr = d_crs_view.row_ptr();
                auto d_colidx = d_crs_view.colidx();
                auto d_values = d_crs_view.values();

                auto o_row_ptr = o_crs_view.row_ptr();
                auto o_colidx = o_crs_view.colidx();
                auto o_values = o_crs_view.values();

                // Count incoming number of indices per process
                std::stringstream ss;
                for (SizeType i = 0; i < nghosts; ++i) {
                    int r = find_rank(ghosts[i]);
                    ++counter[r];
                    // ss << ghosts[i] << " -> " << r << "\n";
                }

                for (int r = 0; r < comm_size; ++r) {
                    if (r == comm.rank() || counter[r] == 0) continue;

                    recv_list[r].resize(counter[r]);
                    // ss << "recevies " << counter[r] << " indices from " << r << "\n";
                }

                // Count outgoing number of indices per process
                std::fill(std::begin(counter), std::end(counter), 0);

                std::vector<SizeType> encountered_ranks;
                std::vector<bool> is_rank_counted(comm.size(), false);

                SizeType local_rows = o_crs_view.rows();
                for (SizeType i = 0; i < n_interface; ++i) {
                    SizeType original_local_idx = local_interface_idx[i];
                    auto row_begin = o_row_ptr[original_local_idx];
                    auto row_end = o_row_ptr[original_local_idx + 1];

                    auto n_cols = row_end - row_begin;

                    if (encountered_ranks.size() < n_cols) {
                        encountered_ranks.resize(n_cols, -1);
                    }

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType col = o_colidx[k];
                        int rank = find_rank(ghosts[col]);

                        if (!is_rank_counted[rank]) {
                            ++counter[rank];
                        }

                        is_rank_counted[rank] = true;
                        encountered_ranks[k - row_begin] = rank;
                    }

                    for (SizeType k = 0; k < n_cols; ++k) {
                        assert(encountered_ranks[k] >= 0);
                        is_rank_counted[encountered_ranks[k]] = false;
                    }

#ifndef NDEBUG
                    for (auto irc : is_rank_counted) {
                        assert(!irc);
                    }
#endif
                }

                for (int r = 0; r < comm_size; ++r) {
                    if (r == comm.rank() || counter[r] == 0) continue;

                    send_list[r].reserve(counter[r]);

                    // ss << "sends " << counter[r] << " indices to " << r << "\n";
                }

                // Fill send_list with global indices

                // ss << "---------------------------\n";
                for (SizeType i = 0; i < n_interface; ++i) {
                    SizeType original_local_idx = local_interface_idx[i];
                    SizeType new_global_idx = interface_offset + i;

                    auto row_begin = o_row_ptr[original_local_idx];
                    auto row_end = o_row_ptr[original_local_idx + 1];
                    auto n_cols = row_end - row_begin;

                    // ss << i << " -> " << original_local_idx << " -> " << new_global_idx << " | " << n_cols << "\n";

                    if (encountered_ranks.size() < n_cols) {
                        encountered_ranks.resize(n_cols, -1);
                    }

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType col = o_colidx[k];
                        int rank = find_rank(ghosts[col]);

                        if (!is_rank_counted[rank]) {
                            send_list[rank].push_back(new_global_idx);
                        }

                        is_rank_counted[rank] = true;
                        encountered_ranks[k - row_begin] = rank;
                    }

                    for (SizeType k = 0; k < n_cols; ++k) {
                        assert(encountered_ranks[k] >= 0);
                        is_rank_counted[encountered_ranks[k]] = false;
                    }

#ifndef NDEBUG
                    for (auto irc : is_rank_counted) {
                        assert(!irc);
                    }
#endif
                }

                // for (int r = 0; r < comm_size; ++r) {
                //     if (send_list[r].empty()) continue;

                //     ss << "to rank " << r << " we send: ";
                //     for (auto idx : send_list[r]) {
                //         ss << idx << " ";
                //     }

                //     ss << "\n";
                // }

                for (int r = 0; r < comm_size; ++r) {
                    if (send_list[r].empty()) {
                        assert(recv_list[r].empty());
                        continue;
                    }

                    int tag = 0;
                    MPI_Status status;
                    auto &send_buff = send_list[r];
                    auto &recv_buff = recv_list[r];

                    assert(!recv_buff.empty());

                    MPI_Sendrecv(send_buff.data(),
                                 send_buff.size(),
                                 MPIType<SizeType>::value(),
                                 r,
                                 tag,
                                 recv_buff.data(),
                                 recv_buff.size(),
                                 MPIType<SizeType>::value(),
                                 r,
                                 tag,
                                 comm.raw_comm(),
                                 &status);
                }

                // for (int r = 0; r < comm_size; ++r) {
                //     if (recv_list[r].empty()) continue;

                //     ss << "from rank " << r << " we received: ";
                //     for (auto idx : recv_list[r]) {
                //         ss << idx << " ";
                //     }

                //     ss << "\n";
                // }

                SizeType linear_index = 0;
                for (int r = 0; r < comm_size; ++r) {
                    if (recv_list[r].empty()) continue;

                    for (auto idx : recv_list[r]) {
                        // ss << linear_index << " -> " << idx << " (old: " << ghosts[linear_index] << ")\n";
                        interface_ghosts[linear_index++] = idx;
                    }
                }

                ////////////////////////////////////////////////////////////////////////////////

                auto G_layout = layout(comm, n_interface, Traits::determine());
                auto GG_layout = square_matrix_layout(G_layout);

                IndexSet d_nnz(n_interface, 0), o_nnz(n_interface, 0);

                {
                    for (SizeType i = 0; i < n_interface; ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = d_row_ptr[original_local_idx];
                        auto row_end = d_row_ptr[original_local_idx + 1];
                        auto n_cols = row_end - row_begin;

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];
                            if (!r.inside(min_idx[col]) || !r.inside(max_idx[col])) {
                                ++d_nnz[i];
                            }
                        }
                    }

                    // count off proc nnz
                    for (SizeType i = 0; i < n_interface; ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = o_row_ptr[original_local_idx];
                        auto row_end = o_row_ptr[original_local_idx + 1];
                        auto n_cols = row_end - row_begin;
                        o_nnz[i] = n_cols;
                    }
                }

                ////////////////////////////////////////////////////////////////////////////////

                Matrix A_GG;
                Vector b_G;
                A_GG.sparse(GG_layout, d_nnz, o_nnz);
                b_G.zeros(G_layout);

                {
                    Write<Matrix> w_A(A_GG);
                    Write<Vector> w_b(b_G);

                    auto b_view = local_view_device(b);

                    for (SizeType i = 0; i < n_interface; ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = d_row_ptr[original_local_idx];
                        auto row_end = d_row_ptr[original_local_idx + 1];
                        auto n_cols = row_end - row_begin;
                        SizeType new_row = interface_offset + i;

                        b_G.set(new_row, b_view.get(original_local_idx));

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];
                            SizeType new_col_local = global_to_local[col];

                            if (new_col_local != -1) {
                                SizeType new_col = interface_offset + new_col_local;
                                Scalar value = d_values[k];
                                A_GG.set(new_row, new_col, value);
                            }
                        }
                    }

                    // count off proc nnz
                    for (SizeType i = 0; i < n_interface; ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = o_row_ptr[original_local_idx];
                        auto row_end = o_row_ptr[original_local_idx + 1];
                        auto n_cols = row_end - row_begin;
                        SizeType new_row = interface_offset + i;

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            Scalar value = o_values[k];
                            SizeType col = o_colidx[k];
                            SizeType new_col = interface_ghosts[col];
                            A_GG.set(new_row, new_col, value);
                        }
                    }
                }

                ////////////////////////////////////////////////////////////////////////////////

                Matrix A_GI, A_IG;
                Vector b_I;

                auto GI_layout = serial_layout(n_interface, n_local);
                auto IG_layout = serial_layout(n_local, n_interface);
                auto I_layout = serial_layout(n_local);

                std::fill(std::begin(d_nnz), std::end(d_nnz), 0);
                std::fill(std::begin(o_nnz), std::end(o_nnz), 0);

                {
                    for (SizeType i = 0; i < n_interface; ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = d_row_ptr[original_local_idx];
                        auto row_end = d_row_ptr[original_local_idx + 1];
                        auto n_cols = row_end - row_begin;

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];
                            if (r.inside(min_idx[col]) && r.inside(max_idx[col])) {
                                ++d_nnz[i];
                            }
                        }
                    }

                    A_GI.sparse(GI_layout, d_nnz, o_nnz);

                    Write<Matrix> w(A_GI);

                    for (SizeType i = 0; i < n_interface; ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = d_row_ptr[original_local_idx];
                        auto row_end = d_row_ptr[original_local_idx + 1];
                        auto n_cols = row_end - row_begin;

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];
                            if (r.inside(min_idx[col]) && r.inside(max_idx[col])) {
                                Scalar value = d_values[k];
                                A_GI.set(i, col, value);
                            }
                        }
                    }
                }

                /////////////////////////////////////////////////

                d_nnz.resize(n_local);
                o_nnz.resize(n_local);

                std::fill(std::begin(d_nnz), std::end(d_nnz), 0);
                std::fill(std::begin(o_nnz), std::end(o_nnz), 0);

                {
                    for (SizeType i = 0; i < n_local; ++i) {
                        auto row_begin = d_row_ptr[i];
                        auto row_end = d_row_ptr[i + 1];
                        auto n_cols = row_end - row_begin;

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];

                            if (!r.inside(min_idx[col]) || !r.inside(max_idx[col])) {
                                ++d_nnz[i];
                            }
                        }
                    }

                    A_IG.sparse(IG_layout, d_nnz, o_nnz);
                    b_I.zeros(I_layout);

                    Write<Matrix> w_A(A_IG);
                    Write<Vector> w_b(b_I);

                    auto b_view = local_view_device(b);

                    for (SizeType i = 0; i < n_local; ++i) {
                        if (global_to_local[i] != -1) continue;

                        auto row_begin = d_row_ptr[i];
                        auto row_end = d_row_ptr[i + 1];
                        auto n_cols = row_end - row_begin;

                        b_I.set(i, b_view.get(i));

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];

                            SizeType new_col_local = global_to_local[col];

                            if (new_col_local != -1) {
                                Scalar value = d_values[k];
                                A_IG.set(i, new_col_local, value);
                            }
                        }
                    }
                }

                ////////////////////////////////////////////////////////////////////////////////

                // comm.synched_print(ss.str());
                // comm.barrier();

                // disp(A_GG);
                // disp(b_G);

                // for (SizeType r = 0; r < comm.size(); ++r) {
                //     comm.barrier();

                //     if (r == comm.rank()) {
                //         std::cout << "=========================================\n";
                //         std::cout << "Rank: " << r << "\n";
                //         std::cout << "=========================================\n";
                //         std::cout << "A_II:\n";
                //         disp(A_II);

                //         disp(b_I);

                //         std::cout << "A_IG:\n";
                //         disp(A_IG);

                //         std::cout << "A_GI:\n";
                //         disp(A_GI);
                //         std::cout << "\n" << std::flush;
                //     }
                // }

                // comm.barrier();

                BDDOperator<Matrix, Vector> op;
                op.init(make_ref(A_GG), make_ref(A_GI), make_ref(A_II), make_ref(A_IG));
                op.init_rhs(b_G, b_I);

                c.stop();
                c_ss << "BDDOperator initialization\n" << c << "\n";

                c.start();

                ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
                cg.verbose(verbose);
                cg.atol(1e-18);
                cg.rtol(1e-18);
                cg.stol(1e-18);

                Vector x_G(row_layout(A_GG), 1);
                cg.solve(op, *op.secant_G_, x_G);

                Vector x_I;
                op.finalize(x_G, b_I, x_I);

                Vector x(layout(b), -666);

                {
                    auto x_view = local_view_device(x);
                    auto x_I_view = local_view_device(x_I);
                    auto x_G_view = local_view_device(x_G);

                    parallel_for(
                        local_range_device(x), UTOPIA_LAMBDA(const SizeType i) { x_view.set(i, x_I_view.get(i)); });

                    parallel_for(local_range_device(x_G), [&](const SizeType i) {
                        SizeType i_local = local_to_global[i] - r.begin();
                        x_view.set(i_local, x_G_view.get(i));
                    });
                }

                c.stop();
                c_ss << "solve\n" << c << "\n";

                if (verbose) {
                    x.comm().root_print(c_ss.str());
                }

                // rename("b", b);
                // write("B.m", b);

                // rename("secant", *op.secant_G_);
                // write("S.m", *op.secant_G_);

                // rename("b_G", b_G);
                // write("B_G.m", b_G);

                // rename("b_I_" + std::to_string(comm.rank()), b_I);
                // write("B_I_" + std::to_string(comm.rank()) + ".m", b_I);

                if (!empty(oracle)) {
                    Scalar diff = norm1(x - oracle);

                    if (diff > 1e-6 && verbose) {
                        comm.root_print(diff);

                        rename("o", oracle);
                        write("O.m", oracle);

                        rename("x", x);
                        write("X.m", x);
                    }

                    utopia_test_assert(diff < 1e-6);
                }
            }

            // comm.barrier();
            // Utopia::Abort("BYE!");
        }

        void poly_qp() {
            QPSolverTest<Matrix, Vector> s;

            OmniQPSolver<Matrix, Vector> solver;
            InputParameters in;

            in.set("backend", "any");
            in.set("type", "pg");
            // in.set("verbose", true);
            in.set("max-it", 2000);

            solver.read(in);

            s.run_qp_solver(solver);
        }
    };

    template <class Matrix, class Vector>
    class MonotoneMGTest {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

        void monotone_amg_test() {
            // Does not work
            const static bool verbose = true;
            const static bool use_masks = false;
            int n_levels = 4;
            int n_coarse = 50;

            auto qp_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            auto fine_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>();
            auto agglomerator = std::make_shared<Agglomerate<Matrix>>();

            MonotoneAlgebraicMultigrid<Matrix, Vector> amg(qp_smoother, fine_smoother, direct_solver, agglomerator);
            amg.set_n_levels(n_levels);

            using ProblemType = utopia::Poisson1D<Matrix, Vector>;
            MultiLevelTestProblem1D<Matrix, Vector, ProblemType> ml_problem(n_levels, n_coarse, !use_masks);
            auto funs = ml_problem.get_functions();

            Vector x, g;
            Matrix H;

            funs.back()->get_eq_constrains_values(x);
            funs.back()->gradient(x, g);
            funs.back()->hessian(x, H);

            Vector lower_bound(layout(g), -0.8), upper_bound(layout(g), 200.);

            amg.verbose(verbose);
            amg.set_box_constraints(make_box_constaints(make_ref(lower_bound), make_ref(upper_bound)));

            utopia_test_assert(amg.solve(H, g, x));

            rename("x", x);
            write("X.m", x);
        }

        void monotone_mg_test() {
            const std::string data_path = Utopia::instance().get("data_path");

            const static bool verbose = false;
            const static bool use_masks = false;

            int n_levels = 6;
            int n_coarse = 50;

            using ProblemType = utopia::Poisson1D<Matrix, Vector>;
            MultiLevelTestProblem1D<Matrix, Vector, ProblemType> ml_problem(n_levels, n_coarse, !use_masks);
            auto funs = ml_problem.get_functions();

            Vector x, g;
            Matrix H;

            funs.back()->get_eq_constrains_values(x);
            funs.back()->gradient(x, g);
            funs.back()->hessian(x, H);

            auto fine_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();

            auto coarse_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            // auto coarse_smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();
            // auto coarse_smoother = std::make_shared<KSPSolver<Matrix, Vector>>();
            // coarse_smoother->pc_type("bjacobi");
            // coarse_smoother->ksp_type("cg");

            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>();
            // auto direct_solver = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();

            MonotoneMultigrid<Matrix, Vector> multigrid(fine_smoother, coarse_smoother, direct_solver);

            std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> interpolation_operators;
            interpolation_operators.resize(n_levels - 1);

            auto &transfers = ml_problem.get_transfer();
            for (SizeType i = 0; i < n_levels - 2; ++i) {
                interpolation_operators[i] = transfers[i];
            }

            auto t = std::static_pointer_cast<MatrixTransfer<Matrix, Vector>>(transfers[n_levels - 2]);
            interpolation_operators[n_levels - 2] = multigrid.new_fine_level_transfer(std::make_shared<Matrix>(t->I()));

            // std::make_shared<IPRTruncatedTransfer<Matrix, Vector>>(std::make_shared<Matrix>(t->I()));

            Vector lower_bound(layout(g), -0.8), upper_bound(layout(g), 200.);

            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.max_it(40);
            multigrid.pre_smoothing_steps(5);
            multigrid.post_smoothing_steps(5);
            multigrid.verbose(verbose);
            multigrid.set_box_constraints(make_box_constaints(make_ref(lower_bound), make_ref(upper_bound)));
            multigrid.update(make_ref(H));

            // avoids flip-floping of active nodes (and Galerkin assembly when nothing changes)
            multigrid.active_set().tol(1e-15);
            multigrid.apply(g, x);

            // disp(x);
        }

        void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(monotone_mg_test);

            // FIXME
            // UTOPIA_RUN_TEST(monotone_amg_test);
        }
    };

    static void qp_solver() {
#ifdef UTOPIA_WITH_PETSC
        PQPSolverTest<PetscMatrix, PetscVector>().run();
        QPSolverTest<PetscMatrix, PetscVector>().run();
        QPSolverTest<PetscMatrix, PetscVector>().run_GS_QR();

        MonotoneMGTest<PetscMatrix, PetscVector>().run();
        // ProjectedGaussSeidelNewTest<PetscMatrix, PetscVector>().run();

#endif  // UTOPIA_WITH_PETSC

#ifdef UTOPIA_WITH_TRILINOS
        QPSolverTest<TpetraMatrixd, TpetraVectord>().run();
#endif  // UTOPIA_WITH_TRILINOS

#ifdef UTOPIA_WITH_BLAS
        QPSolverTest<BlasMatrixd, BlasVectord>().run();  // TODO(zulianp): : because blas is missing min operation ....
#endif                                                   // UTOPIA_WITH_BLAS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(qp_solver);
}  // namespace utopia

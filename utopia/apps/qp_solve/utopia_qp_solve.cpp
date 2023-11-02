#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_PETSC

#include "utopia.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_petsc.hpp"

#include "utopia_MPITimeStatistics.hpp"

#include "utopia_MonotoneMultigrid.hpp"
#include "utopia_MultigridQR.hpp"
#include "utopia_polymorphic_QPSolver.hpp"

#include "utopia_ShiftedPenaltyQPSolver_impl.hpp"

// std
#include <cmath>

namespace utopia {

    template <class Matrix, class Vector>
    class QPSolveApp {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using Communicator = typename Traits<Vector>::Communicator;

        class Problem : public Configurable {
        public:
            void read(Input &in) override {
                in.require("A", A);
                in.require("T", T);

                in.require("A", A);
                in.require("O", O);
                in.require("rhs", rhs);
                in.require("upper_bound", upper_bound);
                in.require("is_contact", is_contact);

                std::stringstream ss;
                ss << "A: " << A.local_size() << " " << A.size() << "\n";
                ss << "T: " << T.local_size() << " " << T.size() << "\n";
                ss << "O: " << O.local_size() << " " << O.size() << "\n";
                ss << "rhs: " << rhs.local_size() << " " << rhs.size() << "\n";
                ss << "is_contact: " << is_contact.local_size() << " " << is_contact.size() << "\n";
                ss << "upper_bound: " << upper_bound.local_size() << " " << upper_bound.size() << "\n";

                A.comm().synched_print(ss.str());
            }

            Matrix A, T, O;
            Vector rhs, upper_bound, is_contact;
        };

        void solve_with_shifted_penalty(Input &in) {
            Problem problem;
            problem.read(in);

            ShiftedPenaltyQPSolver<Matrix> qp_solver;
            qp_solver.read(in);
            auto ptrafo = std::make_shared<PrototypeTransformation<Matrix, Vector>>();

            ptrafo->fun_transform_position = [&](const Vector &x, Vector &x_trafo) {
                Vector x_temp = problem.T * x;
                x_trafo = problem.O * x_temp;
            };

            ptrafo->fun_inverse_transform_direction = [&](const Vector &grad_trafo, Vector &grad) {
                Vector grad_temp;
                grad_temp = problem.O * grad_trafo;
                grad = transpose(problem.T) * grad_temp;
            };

            qp_solver.set_transform(ptrafo);
            BoxConstraints<Vector> box(nullptr, make_ref(problem.upper_bound));
            qp_solver.set_box_constraints(box);

            Vector x(layout(problem.rhs), 0);
            bool ok = qp_solver.solve(problem.A, problem.rhs, x);
        }
    };

    void qp_solve(Input &in) {
        QPSolveApp<PetscMatrix, PetscVector> app;
        app.solve_with_shifted_penalty(in);
    }

    UTOPIA_REGISTER_APP(qp_solve);

    template <class Matrix, class Vector>
    class MGQPSolveApp {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using Communicator = typename Traits<Vector>::Communicator;

        class Problem : public Configurable {
        public:
            void read(Input &in) override {
                const std::string data_path = Utopia::instance().get("data_path") + "/forQR";
                bool use_QR = false;

                Path dir = data_path;
                in.get("dir", dir);
                in.get("verbose", verbose);

                Path rhs_path = dir / "rhs";
                in.get("rhs", rhs_path);
                if (!utopia::read(rhs_path, rhs)) {
                    Utopia::Abort("Could not read r " + rhs_path.to_string() + "!");
                }

                if (verbose) {
                    utopia::out() << "size(rhs): " << size(rhs) << "\n";
                }

                Path x_path = dir / "x";
                in.get("x", x_path);
                if (!utopia::read(x_path, x)) {
                    Utopia::Abort("Could not read " + x_path.to_string() + "!");
                }

                if (verbose) {
                    utopia::out() << "size(x): " << size(x) << "\n";
                }

                Path A_path = dir / "A";
                in.get("A", A_path);
                if (!utopia::read(A_path, A)) {
                    Utopia::Abort("Could not read " + A_path.to_string() + "!");
                }

                if (verbose) {
                    utopia::out() << "size(A): " << size(A) << "\n";
                }

                if (use_QR) {
                    Path Q_path = dir / "Q";
                    in.get("Q", Q_path);

                    if (!Q_path.empty() && utopia::read(Q_path, Q)) {
                        if (verbose) {
                            utopia::out() << "using Q matrix\n";
                        }
                    }

                    if (verbose) {
                        utopia::out() << "size(Q): " << size(Q) << "\n";
                    }

                    Path R_path = dir / "R";
                    in.get("R", R_path);
                    if (!R_path.empty() && utopia::read(R_path, R)) {
                        if (verbose) {
                            utopia::out() << "using R matrix\n";
                        }
                    }

                    if (verbose) {
                        utopia::out() << "size(R): " << size(R) << "\n";
                    }

                    bool transpose_R = true;
                    in.get("transpose_R", transpose_R);

                    if (transpose_R) {
                        R = transpose(R);
                    }
                }

                Path Rot_path = dir / "Rot";
                in.get("Rot", Rot_path);
                if (!utopia::read(Rot_path, Rot)) {
                    Utopia::Abort("Could not read R " + Rot_path.to_string() + "!");
                }

                if (verbose) {
                    utopia::out() << "size(Rot): " << size(Rot) << "\n";
                }

                Path ub_path = dir / "ub";
                in.get("ub", ub_path);
                if (!utopia::read(ub_path, ub)) {
                    Utopia::Abort("Could not read  " + ub_path.to_string() + "!");
                }

                if (verbose) {
                    utopia::out() << "size(ub): " << size(ub) << "\n";
                }

                Path lb_path = dir / "lb";
                in.get("lb", lb_path);
                if (!utopia::read(lb_path, lb)) {
                    Utopia::Abort("Could not read  " + lb_path.to_string() + "!");
                }

                if (verbose) {
                    utopia::out() << "size(lb): " << size(lb) << "\n";
                }

                Path interpolation_base_name = dir / "I";

                in.get("I", interpolation_base_name);

                int n_interpolators = 4;
                int interp_index_offset = 1;
                in.get("n_interpolators", n_interpolators);
                in.get("interp_index_offset", interp_index_offset);

                interpolation_operators.resize(n_interpolators);

                for (int l = 0; l < n_interpolators; ++l) {
                    Path interp_path = interpolation_base_name + std::to_string(l + interp_index_offset);

                    auto mat = std::make_shared<Matrix>();
                    if (!utopia::read(interp_path, *mat)) {
                        Utopia::Abort("Could not read " + interp_path.to_string());
                    }

                    if (verbose) {
                        utopia::out() << "size(I" + std::to_string(l + interp_index_offset) + "): " << size(*mat)
                                      << "\n";
                    }

                    if (l == n_interpolators - 1) {
                        if (use_QR_mg()) {
                            Matrix temp = transpose(Q) * Rot * (*mat);
                            *mat = std::move(temp);
                        } else {
                            Matrix temp = Rot * (*mat);
                            *mat = std::move(temp);
                        }

                        interpolation_operators[l] = std::make_shared<IPTruncatedTransfer<Matrix, Vector>>(mat);
                    } else {
                        interpolation_operators[l] = std::make_shared<IPTransfer<Matrix, Vector>>(mat);
                    }
                }
            }

            void create_linear_system(Matrix &A, Vector &rhs, Vector &x) {
                if (use_QR_mg()) {
                    A = transpose(Q) * Rot * this->A * Rot * Q;
                    rhs = transpose(Q) * Rot * this->rhs;
                    x = transpose(Q) * Rot * this->x;
                } else {
                    A = Rot * this->A * Rot;
                    rhs = Rot * this->rhs;
                    x = Rot * this->x;
                }
            }

            void report_result(Vector &x) {
                if (use_QR_mg()) {
                    this->x = Rot * Q * x;
                } else {
                    this->x = Rot * x;
                }
            }

            bool use_QR_mg() { return !empty(Q) && !empty(R); }

            bool verbose{true};
            Vector rhs, x;
            Vector ub, lb;
            Matrix A, R, Q, Rot;

            std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> interpolation_operators;
        };

        void run_mg_qr(Input &in) {
            Communicator comm;
            MPITimeStatistics stats(comm);

            Problem problem;
            problem.read(in);

            stats.stop_collect_and_restart("init");

            if (problem.use_QR_mg()) {
                solve_problem_QR(in, problem);
            } else {
                solve_problem(in, problem);
            }

            stats.stop_collect_and_restart("solve");

            stats.describe(utopia::out().stream());
        }

        void solve_problem(Input &in, Problem &problem) {
            Matrix QtAQ;
            Vector Qtrhs;
            Vector Qtx;

            problem.create_linear_system(QtAQ, Qtrhs, Qtx);

            auto fine_smoother = std::make_shared<OmniQPSolver<Matrix, Vector>>();
            fine_smoother->set("any", "pgs");
            in.get("qp_smoother", *fine_smoother);

            // auto coarse_smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();
            auto coarse_smoother = fine_smoother;
            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>("mumps", "lu");
            MonotoneMultigrid<Matrix, Vector> multigrid(fine_smoother, coarse_smoother, direct_solver);

            multigrid.set_transfer_operators(problem.interpolation_operators);
            multigrid.max_it(40);
            multigrid.pre_smoothing_steps(3);
            multigrid.post_smoothing_steps(3);
            multigrid.verbose(problem.verbose);

            multigrid.read(in);

            multigrid.set_box_constraints(make_box_constaints(make_ref(problem.lb), make_ref(problem.ub)));

            multigrid.solve(QtAQ, Qtrhs, Qtx);
            problem.report_result(Qtx);
        }

        void solve_problem_QR(Input &in, Problem &problem) {
            Matrix QtAQ;
            Vector Qtrhs;
            Vector Qtx;

            problem.create_linear_system(QtAQ, Qtrhs, Qtx);

            auto fine_smoother = std::make_shared<ProjectedGaussSeidelQR<Matrix, Vector>>();
            fine_smoother->set_R(problem.R);

            auto coarse_smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();
            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>("mumps", "lu");
            MonotoneMultigrid<Matrix, Vector> multigrid(fine_smoother, coarse_smoother, direct_solver);

            multigrid.set_transfer_operators(problem.interpolation_operators);
            multigrid.max_it(40);
            multigrid.pre_smoothing_steps(3);
            multigrid.post_smoothing_steps(3);
            multigrid.verbose(problem.verbose);

            multigrid.read(in);

            multigrid.set_box_constraints(make_box_constaints(make_ref(problem.lb), make_ref(problem.ub)));

            multigrid.solve(QtAQ, Qtrhs, Qtx);
            problem.report_result(Qtx);
        }
    };

    void mg_qp_solve(Input &in) {
        MGQPSolveApp<PetscMatrix, PetscVector> app;
        app.run_mg_qr(in);
    }

    UTOPIA_REGISTER_APP(mg_qp_solve);

}  // namespace utopia

#endif  // UTOPIA_WITH_PETSC

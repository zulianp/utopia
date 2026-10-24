#include "utopia_Testing.hpp"

#include "utopia.hpp"

#include "utopia_LogBarrierFactory.hpp"
#include "utopia_LogBarrierFunction.hpp"
#include "utopia_QuadraticFunction.hpp"
#include "utopia_ShiftedPenalty.hpp"
#include "utopia_ShiftedPenaltyQPSolver_impl.hpp"
#include "utopia_TransformedBoxConstraints.hpp"

#ifdef UTOPIA_ENABLE_PETSC
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"
#endif  // UTOPIA_ENABLE_PETSC

#include <cmath>
#include <string>

namespace utopia {

    // Consistency of value, gradient and hessian of the Penalty implementations (LogBarrier, BoundedLogBarrier,
    // ShiftedPenalty), with and without a transformation to the constraint coordinates
    template <class Matrix, class Vector>
    class PenaltyTest {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Layout = typename Traits::Layout;
        using Comm = typename Traits::Communicator;
        using Penalty = utopia::Penalty<Matrix, Vector>;
        using Transformation = utopia::Transformation<Matrix, Vector>;

        // Even, so that the 2x2 blocks of the transformation are local
        static constexpr SizeType n_local = 12;

        Layout vl;
        std::shared_ptr<Matrix> T;
        std::shared_ptr<Transformation> transformation;

        PenaltyTest() {
            auto &&comm = Comm::get_default();
            vl = layout(comm, n_local, n_local * comm.size());
            init_transformation();
        }

        template <class F>
        void fill(Vector &v, F f) const {
            v.zeros(vl);
            Write<Vector> w(v);
            auto r = range(v);
            for (auto i = r.begin(); i < r.end(); ++i) {
                v.set(i, f(i));
            }
        }

        // Block diagonal with 2x2 rotations and reflections (as the normal-tangential frames of contact)
        void init_transformation() {
            T = std::make_shared<Matrix>();
            T->sparse(square_matrix_layout(vl), 2, 0);

            {
                Write<Matrix> w(*T);
                auto r = row_range(*T);
                for (auto i = r.begin(); i < r.end(); ++i) {
                    const SizeType block = i / 2;
                    const SizeType first = 2 * block;
                    const Scalar theta = 0.3 + 0.7 * block;
                    const Scalar c = std::cos(theta), s = std::sin(theta);
                    const bool reflection = block % 3 == 0;

                    if (i == first) {
                        T->set(i, first, c);
                        T->set(i, first + 1, reflection ? s : -s);
                    } else {
                        T->set(i, first, s);
                        T->set(i, first + 1, reflection ? -c : c);
                    }
                }
            }

            auto t = std::make_shared<PrototypeTransformation<Matrix, Vector>>();
            auto T_ptr = T;
            t->fun_transform_position = [T_ptr](const Vector &in, Vector &out) { out = (*T_ptr) * in; };
            t->fun_transform_direction = [T_ptr](const Vector &in, Vector &out) { out = (*T_ptr) * in; };
            t->fun_inverse_transform_position = [T_ptr](const Vector &in, Vector &out) {
                out = transpose(*T_ptr) * in;
            };
            t->fun_inverse_transform_direction = [T_ptr](const Vector &in, Vector &out) {
                out = transpose(*T_ptr) * in;
            };
            t->fun_transform_matrix = [T_ptr](const Matrix &in, Matrix &out) {
                Matrix temp = in * (*T_ptr);
                out = transpose(*T_ptr) * temp;
            };

            transformation = t;
        }

        struct Case {
            std::string function_type;
            std::string barrier_subtype{"default"};
            bool upper{true};
            bool transform{false};

            std::string name() const {
                return function_type + (barrier_subtype == "default" ? "" : "(" + barrier_subtype + ")") +
                       (upper ? ", upper bound" : ", lower bound") + (transform ? ", transformed" : "");
            }
        };

        // Bounds and a point x whose constraint coordinates y = T x are at distance 0.05..0.35 from the bound,
        // on the feasible side (barriers) or on the infeasible side (penalty)
        void make_problem(const Case &c, BoxConstraints<Vector> &box, Vector &x) const {
            const bool penetrate = c.function_type == "ShiftedPenalty";

            auto bound = std::make_shared<Vector>();
            fill(*bound, [](SizeType i) { return 0.2 + 0.1 * std::sin(1.3 * i); });

            Vector y;
            fill(y, [&](SizeType i) {
                const Scalar b = 0.2 + 0.1 * std::sin(1.3 * i);
                const Scalar d = 0.05 + 0.3 * (0.5 + 0.5 * std::sin(2.1 * i + 0.4));
                const Scalar feasible_direction = c.upper ? -1 : 1;
                return b + (penetrate ? -d : d) * feasible_direction;
            });

            if (c.upper) {
                box = BoxConstraints<Vector>(nullptr, bound);
            } else {
                box = BoxConstraints<Vector>(bound, nullptr);
            }

            if (c.transform) {
                x = transpose(*T) * y;
            } else {
                x = y;
            }
        }

        std::unique_ptr<Penalty> make_penalty(const Case &c, BoxConstraints<Vector> &box) const {
            auto p = PenaltyFactory<Matrix, Vector>::new_penalty(c.function_type);

            InputParameters params;
            params.set("barrier_parameter", 1.);
            params.set("min_barrier_parameter", 1.);
            params.set("barrier_thickness", 0.5);
            params.set("barrier_subtype", c.barrier_subtype);
            params.set("penalty_parameter", 10.);
            p->read(params);

            p->set_box_constraints(make_ref(box));
            if (c.transform) {
                p->set_transform(transformation);
            }

            return p;
        }

        Scalar value(const Penalty &p, const Vector &x) const {
            Scalar v = 0;
            utopia_test_assert(p.value(x, v));
            return v;
        }

        Vector gradient(const Penalty &p, const Vector &x) const {
            Vector g(vl, 0.);
            utopia_test_assert(p.gradient(x, g));
            return g;
        }

        void check(const Case &c) const {
            BoxConstraints<Vector> box;
            Vector x;
            make_problem(c, box, x);
            auto p = make_penalty(c, box);

            // Barrier continuation steps (no-op for min_barrier_parameter = barrier_parameter) or shift updates.
            // After two updates the shift of ShiftedPenalty differs from the violation
            p->update(x);
            p->update(x);

            Vector d;
            fill(d, [](SizeType i) { return std::cos(0.7 * i + 0.2); });

            const Scalar eps = 1e-6;
            Vector x_plus = x + eps * d;
            Vector x_minus = x - eps * d;

            // value() adds to the given value
            Scalar v_offset = 1;
            utopia_test_assert(p->value(x, v_offset));
            const Scalar v = value(*p, x);
            const Scalar accumulate_error = std::abs((v_offset - 1) - v);

            // Directional derivative of the value
            const Vector g = gradient(*p, x);
            const Scalar dv_fd = (value(*p, x_plus) - value(*p, x_minus)) / (2 * eps);
            const Scalar dv = dot(g, d);
            const Scalar gradient_error = std::abs(dv_fd - dv) / std::max(std::abs(dv), Scalar(1e-12));

            // Directional derivative of the gradient
            Matrix H;
            H.identity(square_matrix_layout(vl), 0.);
            utopia_test_assert(p->hessian(x, H));
            Vector Hd = H * d;
            Vector dg_fd = gradient(*p, x_plus) - gradient(*p, x_minus);
            dg_fd *= 1 / (2 * eps);
            const Scalar hessian_error = Scalar(norm_infty(Vector(dg_fd - Hd))) / Scalar(norm_infty(Hd));

            // Diagonal version of the hessian
            Vector h_diag(vl, 0.);
            utopia_test_assert(p->hessian(x, h_diag));
            Vector H_diag = diag(H);
            const Scalar diag_error = Scalar(norm_infty(Vector(h_diag - H_diag))) / Scalar(norm_infty(H_diag));

            if (!x.comm().rank()) {
                utopia::out() << "penalty_functions: " << c.name() << ": |value accumulation| = " << accumulate_error
                              << ", value/gradient = " << gradient_error << ", gradient/hessian = " << hessian_error
                              << ", diag(hessian) = " << diag_error << "\n";
            }

            utopia_test_assert(accumulate_error <= 1e-10 * std::max(Scalar(1), std::abs(v)));
            utopia_test_assert(gradient_error <= 1e-6);
            utopia_test_assert(hessian_error <= 1e-6);
            utopia_test_assert(diag_error <= 1e-12);
        }

        void derivatives() {
            std::vector<Case> cases;
            for (bool transform : {false, true}) {
                for (bool upper : {true, false}) {
                    cases.push_back({"LogBarrier", "default", upper, transform});

                    for (auto subtype : {"default", "polynomial", "polynomial4", "composite_polynomial"}) {
                        cases.push_back({"BoundedLogBarrier", subtype, upper, transform});
                    }

                    cases.push_back({"ShiftedPenalty", "default", upper, transform});
                }
            }

            for (auto &c : cases) {
                check(c);
            }
        }

        // LogBarrierFunction::value is the value of the unconstrained function plus the barrier
        void log_barrier_function_value() {
            Case c{"LogBarrier", "default", true, false};
            BoxConstraints<Vector> box;
            Vector x;
            make_problem(c, box, x);

            auto A = std::make_shared<Matrix>();
            A->identity(square_matrix_layout(vl), 2.);
            auto b = std::make_shared<Vector>(vl, 1.);
            auto fun = std::make_shared<QuadraticFunction<Matrix, Vector>>(A, b);

            auto barrier = std::make_shared<LogBarrier<Matrix, Vector>>(make_ref(box));
            barrier->set_all_barrier_parameters_to(1.);
            LogBarrierFunction<Matrix, Vector> barrier_function(fun, barrier);

            Scalar f = 0, B = 0, F = 0;
            utopia_test_assert(fun->value(x, f));
            utopia_test_assert(barrier->value(x, B));
            utopia_test_assert(barrier_function.value(x, F));

            const Scalar error = std::abs(F - (f + B));

            if (!x.comm().rank()) {
                utopia::out() << "penalty_functions: LogBarrierFunction: f = " << f << ", barrier = " << B
                              << ", f + barrier = " << F << ", error = " << error << "\n";
            }

            utopia_test_assert(error <= 1e-12 * std::abs(f + B));
        }

        // ShiftedPenaltyQPSolver (shift update and Newton step at every iteration) on min 1/2 x^T A x - b^T x,
        // x <= ub, with A = 2 I and half of the constraints active: the solution is min(b / 2, ub)
        void shifted_penalty_qp() {
            Matrix A;
            A.identity(square_matrix_layout(vl), 2.);

            Vector b, ub, expected;
            fill(b, [](SizeType i) { return 1. + 0.5 * std::sin(0.9 * i); });
            fill(ub, [](SizeType i) { return (i % 2 == 0) ? 0.25 : 2.; });
            fill(expected, [](SizeType i) {
                const Scalar b_i = 1. + 0.5 * std::sin(0.9 * i);
                return std::min(b_i / 2, (i % 2 == 0) ? 0.25 : 2.);
            });

            ShiftedPenaltyQPSolver<Matrix, Vector> qp_solver;
            InputParameters params;
            params.set("penalty_parameter", 1e4);
            params.set("max_it", 100);
            params.set("atol", 1e-14);
            params.set("stol", 1e-14);
            qp_solver.read(params);
            qp_solver.set_box_constraints(make_upper_bound_constraints(make_ref(ub)));

            Vector x(vl, 0.);
            const bool converged = qp_solver.solve(A, b, x);
            const Scalar error = norm_infty(Vector(x - expected));

            if (!x.comm().rank()) {
                utopia::out() << "penalty_functions: ShiftedPenaltyQPSolver: converged = " << converged
                              << ", |x - x_exact| = " << error << "\n";
            }

            utopia_test_assert(converged);
            utopia_test_assert(error <= 1e-10);
        }

        void run() {
            UTOPIA_RUN_TEST(derivatives);
            UTOPIA_RUN_TEST(log_barrier_function_value);
            UTOPIA_RUN_TEST(shifted_penalty_qp);
        }
    };

    static void penalty_functions() {
#ifdef UTOPIA_ENABLE_PETSC
        PenaltyTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_ENABLE_PETSC
    }

    UTOPIA_REGISTER_TEST_FUNCTION(penalty_functions);
}  // namespace utopia

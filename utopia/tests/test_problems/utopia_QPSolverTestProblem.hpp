#ifndef UTOPIA_QP_SOLVER_TEST_PROBLEM_HPP
#define UTOPIA_QP_SOLVER_TEST_PROBLEM_HPP

#include <iostream>
#include "utopia_Chrono.hpp"
#include "utopia_Range.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class QPSolverTestProblem {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        template <class QPSolver>
        static void run(const SizeType &n, const bool verbose, QPSolver &qp_solver, const bool use_constraints = true) {
            run(Comm::get_default(), n, verbose, qp_solver, use_constraints);
        }

        template <class QPSolver>
        static void run(const Comm &comm,
                        const SizeType &n,
                        const bool verbose,
                        QPSolver &qp_solver,
                        const bool use_constraints = true) {
            Matrix m;
            m.sparse(layout(comm, Traits::decide(), Traits::decide(), n, n), 3, 2);

            if (verbose) {
                std::cout << comm.rank() << "/" << comm.size() << ": problem size = " << m.local_rows() << std::endl;
            }

            assemble_laplacian_1D(m);
            {
                Range r = row_range(m);
                const SizeType r_begin = r.begin();
                const SizeType r_end = r.end();

                Write<Matrix> w(m);
                if (r_begin == 0) {
                    m.set(0, 0, 1.);
                    m.set(0, 1, 0);
                }

                if (r_end == n) {
                    m.set(n - 1, n - 1, 1.);
                    m.set(n - 1, n - 2, 0);
                }
            }

            Vector rhs(row_layout(m), 1.);
            {
                // Creating test vector (alternative way see [assemble vector alternative], which might be easier for
                // beginners)
                Range r = range(rhs);
                const SizeType r_begin = r.begin();
                const SizeType r_end = r.end();

                Write<Vector> w(rhs);

                if (r_begin == 0) {
                    rhs.set(0, 0);
                }

                if (r_end == n) {
                    rhs.set(n - 1, 0.);
                }
            }

            Vector upper_bound(row_layout(m), 100.0);
            Vector solution(row_layout(m), 0.0);

            utopia_test_assert(upper_bound.comm().size() == rhs.comm().size());

            qp_solver.max_it(n * 40);
            qp_solver.verbose(verbose);

            if (use_constraints) {
                qp_solver.set_box_constraints(make_upper_bound_constraints(make_ref(upper_bound)));
                utopia_test_assert(qp_solver.has_bound());
                utopia_test_assert(qp_solver.get_upper_bound().comm().size() == rhs.comm().size());
            }

            // std::cout << "Solving..." << std::endl;

            Chrono c;
            c.start();
            bool ok = qp_solver.solve(m, rhs, solution);
            c.stop();

            utopia_test_assert(ok);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_QP_SOLVER_TEST_PROBLEM_HPP

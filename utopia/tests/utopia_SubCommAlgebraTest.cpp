#include "test_problems/utopia_QPSolverTestProblem.hpp"
#include "utopia.hpp"
#include "utopia_BlockQPSolver.hpp"
#include "utopia_ParallelTestRunner.hpp"
#include "utopia_SemismoothNewton_old.hpp"
#include "utopia_Testing.hpp"
#include "utopia_assemble_laplacian_1D.hpp"
#include "utopia_petsc_Redundant.hpp"
#include "utopia_petsc_RedundantQPSolver.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class SubCommAlgebraTest final : public AlgebraUnitTest<Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using IndexSet = typename Traits::IndexSet;
        using Comm = typename Traits::Communicator;
        using Layout = typename Traits::Layout;
        using MatrixLayout = typename Traits::MatrixLayout;

        void sum_vectors() {
            Vector v1(layout(this->comm(), 2, this->comm().size() * 2), 1.0);
            Vector v2(layout(v1), 2.0);
            Vector v3 = v1 + v2;
            Scalar norm_v3 = norm1(v3);
            utopia_test_assert(approxeq(norm_v3, this->comm().size() * 3 * 2, device::epsilon<Scalar>()));
        }

        void mat_vec_mult() {
            Vector x(layout(this->comm(), 2, this->comm().size() * 2), 1.0);
            Vector b;

            Matrix A;
            A.sparse(square_matrix_layout(layout(x)), 3, 2);
            assemble_laplacian_1D(A);

            b = A * x;
            utopia_test_assert(b.size() == this->comm().size() * 2);
        }

        void redundant_test() {
            int sub_comms = this->comm().size() > 1 ? 2 : 1;

            if (sub_comms == 1) return;

            Matrix A, A_sub;
            Vector v1, v1_sub, v2_sub;

            auto v_lo = layout(this->comm(), 2, this->comm().size() * 2);
            auto m_lo = square_matrix_layout(v_lo);

            A.sparse(m_lo, 3, 2);
            assemble_laplacian_1D(A);

            v1.zeros(v_lo);
            each_write(v1, [](const SizeType &i) -> Scalar { return Scalar(i); });

            // creating subcommunicators, indices, and scatters
            Redundant<Matrix, Vector> red;
            red.init(v_lo, sub_comms);

            red.create_sub_vector(v1_sub);
            red.create_sub_matrix(A, A_sub);
            red.super_to_sub(v1, v1_sub);

            v2_sub = A_sub * v1_sub;
            const Scalar sub_norm = norm2(v2_sub);

            red.sub_to_super(v2_sub, v1);
            const Scalar super_norm = norm2(v1);

            utopia_test_assert(approxeq(sub_norm, super_norm));
        }

        void solve_problem() {
            {
                MPGRP<Matrix, Vector> solver;
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, solver, true);
            }

            {
                SemismoothNewton<Matrix, Vector> solver(std::make_shared<Factorization<Matrix, Vector>>());
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, solver, true);
            }

            {
                ProjectedGaussSeidel<Matrix, Vector> solver;
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, solver, true);
            }

            {
                SemismoothNewton<Matrix, Vector> solver(std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>());
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, solver, true);
            }

            {
                SemismoothNewton<Matrix, Vector> solver(std::make_shared<MPGRP<Matrix, Vector>>());
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, solver, true);
            }

            // {
            //     SemismoothNewton<Matrix, Vector, INVALID_BACKEND> solver(std::make_shared<MPGRP<Matrix, Vector>>());
            //     QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, true, solver, true);
            // }
        }

        void qp_solver_with_clone() {
            {
                MPGRP<Matrix, Vector> temp_qp;
                auto bqp_ptr = std::shared_ptr<MPGRP<Matrix, Vector>>(temp_qp.clone());
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, *bqp_ptr, true);
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, *bqp_ptr, false);
            }

            {
                ProjectedGaussSeidel<Matrix, Vector> temp_qp;
                temp_qp.use_sweeper(this->comm().size() > 1);
                auto bqp_ptr = std::shared_ptr<ProjectedGaussSeidel<Matrix, Vector>>(temp_qp.clone());
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, *bqp_ptr, true);
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, *bqp_ptr, false);
            }

            {
                auto qp = std::make_shared<SemismoothNewton<Matrix, Vector>>(std::make_shared<MPGRP<Matrix, Vector>>());
                auto bqp_ptr = std::shared_ptr<SemismoothNewton<Matrix, Vector>>(qp->clone());
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, *bqp_ptr, true);
            }

            {
                auto qp = std::make_shared<SemismoothNewton<Matrix, Vector>>(
                    std::make_shared<Factorization<Matrix, Vector>>());
                auto bqp_ptr = std::shared_ptr<SemismoothNewton<Matrix, Vector>>(qp->clone());
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, *bqp_ptr, true);
            }
        }

        void block_qp_solver() {
            {
                auto qp = std::make_shared<MPGRP<Matrix, Vector>>();

                BlockQPSolver<Matrix, Vector> bqp(qp);
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, bqp, true);
            }

            // {
            //     auto qp = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();

            //     BlockQPSolver<Matrix, Vector> bqp(qp);
            //     QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, bqp, true);
            // }

            {
                auto qp = std::make_shared<SemismoothNewton<Matrix, Vector>>(std::make_shared<MPGRP<Matrix, Vector>>());
                // qp->verbose(true);
                BlockQPSolver<Matrix, Vector> bqp(qp);
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, bqp, true);
            }

            {
                auto qp = std::make_shared<SemismoothNewton<Matrix, Vector>>(
                    std::make_shared<Factorization<Matrix, Vector>>());
                // qp->verbose(true);
                BlockQPSolver<Matrix, Vector> bqp(qp);
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, bqp, true);
            }
        }

        void block_qp_solver_with_clone() {
            {
                auto qp = std::make_shared<MPGRP<Matrix, Vector>>();
                BlockQPSolver<Matrix, Vector> temp_qp(qp);
                auto bqp_ptr = std::shared_ptr<BlockQPSolver<Matrix, Vector>>(temp_qp.clone());
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, *bqp_ptr, true);
            }

            // {
            //     auto qp = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            //     qp->use_sweeper(false);

            //     BlockQPSolver<Matrix, Vector> temp_qp(qp);
            //     auto bqp_ptr = std::shared_ptr<BlockQPSolver<Matrix, Vector>>(temp_qp.clone());
            //     QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, true, *bqp_ptr, true);
            // }

            // {
            //     // auto qp = std::make_shared<SemismoothNewton<Matrix, Vector>>(
            //     //     std::make_shared<Factorization<Matrix, Vector>>());

            //     auto qp = std::make_shared<SemismoothNewton<Matrix, Vector>>(std::make_shared<MPGRP<Matrix,
            //     Vector>>());

            //     BlockQPSolver<Matrix, Vector> temp_qp(qp);
            //     auto bqp_ptr = std::shared_ptr<BlockQPSolver<Matrix, Vector>>(temp_qp.clone());
            //     QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, *bqp_ptr, true);
            // }

            {
                auto qp = std::make_shared<SemismoothNewton<Matrix, Vector>>(std::make_shared<MPGRP<Matrix, Vector>>());

                BlockQPSolver<Matrix, Vector> temp_qp(qp);
                auto bqp_ptr = std::shared_ptr<BlockQPSolver<Matrix, Vector>>(temp_qp.clone());
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, *bqp_ptr, true);
            }

            {
                auto qp = std::make_shared<SemismoothNewton<Matrix, Vector>>(
                    std::make_shared<Factorization<Matrix, Vector>>());

                BlockQPSolver<Matrix, Vector> temp_qp(qp);
                auto bqp_ptr = std::shared_ptr<BlockQPSolver<Matrix, Vector>>(temp_qp.clone());
                QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, *bqp_ptr, true);
            }
        }

        void redundant_qp_solver() {
            auto qp = std::make_shared<MPGRP<Matrix, Vector>>();
            // qp->verbose(true);
            RedundantQPSolver<Matrix, Vector> redqp(qp, 2);

            QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, redqp, true);

            // second time should reuse all pre-allocated data
            QPSolverTestProblem<Matrix, Vector>::run(this->comm(), 10, false, redqp, true);
        }

        void run() {
            UTOPIA_RUN_TEST(sum_vectors);
            UTOPIA_RUN_TEST(mat_vec_mult);
            UTOPIA_RUN_TEST(solve_problem);
            UTOPIA_RUN_TEST(redundant_test);
            UTOPIA_RUN_TEST(redundant_qp_solver);
            UTOPIA_RUN_TEST(qp_solver_with_clone);
            UTOPIA_RUN_TEST(block_qp_solver);
            UTOPIA_RUN_TEST(block_qp_solver_with_clone);
        }
    };

    void sub_comm_algebra() {
#ifdef WITH_PETSC
        run_parallel_test<SubCommAlgebraTest<PetscMatrix, PetscVector>>();
#endif  // WITH_PETSC

        // #ifdef WITH_TRILINOS
        //         run_parallel_test< SubCommAlgebraTest<TpetraMatrix, TpetraVector> >();
        // #endif //WITH_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(sub_comm_algebra);

}  // namespace utopia
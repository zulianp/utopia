#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_assemble_laplacian_1D.hpp"
#include "utopia_ParallelTestRunner.hpp"
#include "test_problems/utopia_QPSolverTestProblem.hpp"
#include "utopia_petsc_Redundant.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class SubCommAlgebraTest final : public AlgebraUnitTest<Vector> {
    public:
        using Traits       = utopia::Traits<Vector>;
        using Scalar       = typename Traits::Scalar;
        using SizeType     = typename Traits::SizeType;
        using IndexSet     = typename Traits::IndexSet;
        using Comm         = typename Traits::Communicator;
        using Layout       = typename Traits::Layout;
        using MatrixLayout = typename Traits::MatrixLayout;

        void sum_vectors()
        {
            Vector v1(layout(this->comm(), 2, this->comm().size() * 2), 1.0);
            Vector v2(layout(v1), 2.0);
            Vector v3 = v1 + v2;
            Scalar norm_v3 = norm1(v3);
            utopia_test_assert(approxeq(norm_v3, this->comm().size() * 3 * 2, device::epsilon<Scalar>()));
        }

        void redundant_test()
        {
            int sub_comms = this->comm().size() > 1 ? 2 : 1;

            if(sub_comms == 1) return;

            auto v_lo = layout(this->comm(), 2, this->comm().size() * 2);
            auto m_lo = square_matrix_layout(v_lo);

            Matrix A_sub, A; A.sparse(m_lo, 3, 2);
            assemble_laplacian_1D(A);

            Vector v1(v_lo, 1.0), v1_sub;

            each_write(v1, [](const SizeType &i) -> Scalar {
                return Scalar(i);
            });

            //FIXME
            Redundant<Matrix, Vector> red;
            red.init(v_lo, sub_comms);

            PetscVecScatter sup_to_sub, sub_to_sup;

            red.create_sub_vector(
                        v1,
                        v1_sub,
                        sup_to_sub,
                        sub_to_sup
                    );

            red.create_sub_matrix(A, A_sub);

            disp(A);
            disp(A_sub);

            red.super_to_sub(v1, sup_to_sub, v1_sub);

            // red.super_to_sub(A, A_sub);

            disp(v1);
            disp(v1_sub);

            Vector v2_sub = A_sub * v1_sub;

            // each_write(v1_sub, [](const SizeType &i) -> Scalar {
            //     return Scalar(i);
            // });

            red.super_to_sub(v2_sub, sub_to_sup, v1);

            disp(v1);
            disp(v2_sub);

        }

        void solve_problem()
        {
            MPGRP<Matrix, Vector> solver;
            QPSolverTestProblem<Matrix, Vector>::run(
                this->comm(),
                10,
                false,
                solver,
                true
            );
        }

        void run()
        {
            UTOPIA_RUN_TEST(sum_vectors);
            UTOPIA_RUN_TEST(solve_problem);
            // UTOPIA_RUN_TEST(redundant_test);
        }

    };

    void sub_comm_algebra()
    {
#ifdef WITH_PETSC
        run_parallel_test< SubCommAlgebraTest<PetscMatrix, PetscVector> >();
#endif //WITH_PETSC

// #ifdef WITH_TRILINOS
//         run_parallel_test< SubCommAlgebraTest<TpetraMatrix, TpetraVector> >();
// #endif //WITH_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(sub_comm_algebra);

}
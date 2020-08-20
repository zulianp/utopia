#include "utopia.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_Assert.hpp"
#include "utopia_ParallelTestRunner.hpp"
#include "utopia_Testing.hpp"

#include "utopia_assemble_laplacian_1D.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class Loops final : public AlgebraUnitTest<Vector> {
    private:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using IndexSet = typename Traits::IndexSet;
        using Comm = typename Traits::Communicator;

        // void vector_loop() {}

        void matrix_loop() {
            SizeType n = 5;

            Matrix mat;
            mat.sparse(layout(this->comm(), Traits::decide(), Traits::decide(), n, n), 3, 3);

            // FIXME this is serial and not really portable
            // Scalar norm1_read = 0.0;
            {
                // n x n matrix with maximum 3 entries x row
                Write<Matrix> w(mat);
                Range r = row_range(mat);
                SizeType n = size(mat).get(0);

                const SizeType r_begin = r.begin();
                const SizeType r_end = r.end();

                for (SizeType i = r_begin; i != r_end; ++i) {
                    mat.c_set(i, i, i);
                    mat.c_set(i, n - 1, i);
                    mat.c_set(i, 0, i);
                }
            }

            // mat.read([&](const SizeType &, const SizeType &, const Scalar &val) { norm1_read += device::abs(val); });
            // norm1_read = this->comm().sum(norm1_read);

            const Scalar norm1_mat = sum(abs(mat));
            utopia_test_asserteq(Scalar(26.0), norm1_mat, device::epsilon<Scalar>());
        }

    public:
        void run() override {
            // UTOPIA_RUN_TEST(vector_loop);
            UTOPIA_RUN_TEST(matrix_loop);
        }
    };

    static void loops() {
#ifdef WITH_BLAS
        // Serial backend
        run_serial_test<Loops<BlasMatrixd, BlasVectord>>();
#endif  // WITH_BLAS

#ifdef UTOPIA_WITH_PETSC
        run_parallel_test<Loops<PetscMatrix, PetscVector>>();
#endif  // UTOPIA_WITH_PETSC

#ifdef WITH_TRILINOS
        run_parallel_test<Loops<TpetraMatrixd, TpetraVectord>>();
#endif  // WITH_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(loops);
}  // namespace utopia

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

        void vector_loop() {}

        void matrix_loop() {
            SizeType n = 5;

            Matrix mat;
            mat.sparse(layout(this->comm(), Traits::decide(), Traits::decide(), n, n), 3, 3);

            Scalar norm1_read = 0.0;
            {
                // n x n matrix with maximum 3 entries x row
                Write<Matrix> w(mat);
                Range r = row_range(mat);
                auto n = size(mat).get(0);

                for (SizeType i = r.begin(); i != r.end(); ++i) {
                    mat.c_set(i, i, i);
                    mat.c_set(i, n - 1, i);
                    mat.c_set(i, 0, i);
                }
            }

            std::stringstream ss;
            mat.read([&](const SizeType &i, const SizeType &j, const Scalar &val) { norm1_read += std::abs(val); });

            // this->comm().synched_print(ss.str(), std::cout);

            norm1_read = this->comm().sum(norm1_read);

            const Scalar norm1_mat = sum(abs(mat));

            utopia_test_asserteq(norm1_read, norm1_mat, device::epsilon<Scalar>());
        }

    public:
        void run() override {
            UTOPIA_RUN_TEST(vector_loop);
            UTOPIA_RUN_TEST(matrix_loop);
        }
    };

    static void loops() {
#ifdef WITH_BLAS
        Loops<BlasMatrixd, BlasVectord>().run();
#endif  // WITH_BLAS

#ifdef WITH_PETSC
        Loops<PetscMatrix, PetscVector>().run();
#endif  // WITH_PETSC

#ifdef WITH_TRILINOS
        Loops<TpetraMatrixd, TpetraVectord>().run();
#endif  // WITH_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(loops);
}  // namespace utopia

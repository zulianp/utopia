#include "utopia.hpp"
#include "utopia_Testing.hpp"

#include "utopia_DeprecatedHeaders.hpp"

namespace utopia {

    // All tests deprecated functionalities shall be moved here.
    template <class Matrix, class Vector>
    class DeprecatedTest {
    public:
        using Traits = utopia::Traits<Vector>;
        using SizeType = typename Traits::SizeType;
        using Scalar = typename Traits::Scalar;
        using Comm = typename Traits::Communicator;

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                std::cout << "\nBackend: " << Traits::backend_info().get_name() << std::endl;
            }
        }

        static void vector_factory() {
            Comm comm;

            SizeType n_local = 2;
            SizeType n = n_local * comm.size();

            Vector v;
            v = values(n, 1.0);
            utopia_test_assert(v.size() == n);

            v = local_values(n_local, 1.0);
            utopia_test_assert(v.size() == n);

            v = zeros(n);
            utopia_test_assert(v.size() == n);

            v = local_zeros(n_local);
            utopia_test_assert(v.size() == n);
        }

        static void matrix_factory() {
            Comm comm;

            SizeType n_local = 2;
            SizeType n = n_local * comm.size();

            Matrix m;
            m = sparse(n, n, 3);
            utopia_test_assert(m.rows() == n);
            utopia_test_assert(m.cols() == n);

            m = local_sparse(n_local, n_local, 3);
            utopia_test_assert(m.rows() == n);
            utopia_test_assert(m.cols() == n);

            m = identity(n, n);
            utopia_test_assert(m.rows() == n);
            utopia_test_assert(m.cols() == n);

            m = local_identity(n_local, n_local);
            utopia_test_assert(m.rows() == n);
            utopia_test_assert(m.cols() == n);
        }

        static void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(vector_factory);
            UTOPIA_RUN_TEST(matrix_factory);
        }
    };

    template <class Matrix, class Vector>
    class DenseDeprecatedTest {
    public:
        using Traits = utopia::Traits<Vector>;
        using SizeType = typename Traits::SizeType;
        using Scalar = typename Traits::Scalar;
        using Comm = typename Traits::Communicator;

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                std::cout << "\nBackend: " << Traits::backend_info().get_name() << std::endl;
            }
        }

        static void matrix_factory() {
            Comm comm;

            SizeType n_local = 2;
            SizeType n = n_local * comm.size();

            Matrix m;
            m = dense(n, n);
            utopia_test_assert(m.rows() == n);
            utopia_test_assert(m.cols() == n);

            // m = dense(n_local, n_local);
            // utopia_test_assert(m.rows() == n);
            // utopia_test_assert(m.cols() == n);

            // m = dense_identity(n, n);
            // utopia_test_assert(m.rows() == n);
            // utopia_test_assert(m.cols() == n);

            // m = dense_local_identity(n_local, n_local);
            // utopia_test_assert(m.rows() == n);
            // utopia_test_assert(m.cols() == n);
        }

        static void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(matrix_factory);
        }
    };

    static void deprecated() {
#ifdef WITH_BLAS
        DeprecatedTest<BlasMatrixd, BlasVectord>().run();
        // DenseDeprecatedTest<BlasMatrixd, BlasVectord>().run();
#endif  // WITH_BLAS

#ifdef WITH_PETSC
        DeprecatedTest<utopia::PetscMatrix, utopia::PetscVector>::run();
        // DenseDeprecatedTest<PetscMatrix, PetscVector>().run();
#endif  // WITH_PETSC

#ifdef WITH_TRILINOS
        DeprecatedTest<utopia::TpetraMatrix, utopia::TpetraVector>::run();
#endif  // WITH_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(deprecated);

}  // namespace utopia
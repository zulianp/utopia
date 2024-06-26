#include "utopia_Base.hpp"

#ifdef UTOPIA_ENABLE_DEPRECATED_API

#include "utopia.hpp"
#include "utopia_Testing.hpp"

#include "utopia_DeprecatedHeaders.hpp"

// Remove warnings of depreacted functions since we are using them right now
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

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
                utopia::out() << "\nBackend: " << Traits::backend_info().get_name() << std::endl;
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

            Vector v1 = local_values(n_local, 1.0);
            utopia_test_assert(v1.size() == n);

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
                utopia::out() << "\nBackend: " << Traits::backend_info().get_name() << std::endl;
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
#ifdef UTOPIA_ENABLE_BLAS
        DeprecatedTest<BlasMatrixd, BlasVectord>().run();
        DenseDeprecatedTest<BlasMatrixd, BlasVectord>().run();
#endif  // UTOPIA_ENABLE_BLAS

#ifdef UTOPIA_ENABLE_PETSC
        DeprecatedTest<utopia::PetscMatrix, utopia::PetscVector>::run();
        // DenseDeprecatedTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_ENABLE_PETSC

#ifdef UTOPIA_ENABLE_TRILINOS
        DeprecatedTest<utopia::TpetraMatrix, utopia::TpetraVector>::run();
#endif  // UTOPIA_ENABLE_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(deprecated);

// reset warning to previous stack state
#pragma GCC diagnostic pop

}  // namespace utopia

#endif  // UTOPIA_ENABLE_DEPRECATED_API

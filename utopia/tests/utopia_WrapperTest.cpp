
#include "utopia.hpp"
#include "utopia_Testing.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class WrapperTest {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using IndexSet = typename Traits::IndexSet;
        using Comm = typename Traits::Communicator;

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

        void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(vector_factory_test);
            UTOPIA_RUN_TEST(matrix_factory_test);
            UTOPIA_RUN_TEST(matrix_assembly_test);
        }

        void matrix_assembly_test() {
            auto &&comm = Comm::get_default();

            Matrix mat;
            mat.dense(layout(comm, Traits::decide(), Traits::decide(), n_dofs, n_dofs), 0);

            // Assemble 1D laplacian with dirichlet nodes at the boundary
            {
                Write<Matrix> write(mat);
                Range rr = row_range(mat);

                // FIXME assumed row is owned by this proc completely.
                int rbegin = rr.begin();
                int rend = rr.end();

                if (rbegin == 0) {
                    mat.set(0, 0, 1);
                    rbegin += 1;
                }

                if (rend == mat.size().get(0)) {
                    mat.set(n_dofs - 1, n_dofs - 1, 1);
                    rend -= 1;
                }

                for (int i = rbegin; i != rend; ++i) {
                    mat.set(i, i, 2);
                    mat.set(i, i - 1, -1);
                    mat.set(i, i + 1, -1);
                }
            }

            Vector expected(row_layout(mat), 0.);

            // Assemble expected result
            {
                Write<Vector> write(expected);
                Range r = range(expected);
                if (r.begin() == 0) {
                    expected.set(0, 10);
                }

                if (SizeType(r.end()) == SizeType(expected.size())) {
                    expected.set(expected.size() - 1, 10);
                }
            }

            Vector vec(layout(expected), 10);
            utopia_test_assert(approxeq(expected, mat * vec));
        }

        void matrix_factory_test() {
            auto &&comm = Comm::get_default();
            Matrix mat;
            mat.dense(layout(comm, Traits::decide(), Traits::decide(), n_dofs, n_dofs), 0.1);
            {
                Read<Matrix> read(mat);
                Range rr = row_range(mat);

                for (int i = rr.begin(); i != rr.end(); ++i) {
                    for (int j = 0; j != n_dofs; ++j) {
                        utopia_test_assert(approxeq(mat.get(i, j), 0.1));
                    }
                }
            }

            mat.identity();
            {
                Read<Matrix> read(mat);
                Range rr = row_range(mat);

                for (int i = rr.begin(); i != rr.end(); ++i) {
                    for (int j = 0; j != n_dofs; ++j) {
                        utopia_test_assert(approxeq(mat.get(i, j), i == j));
                    }
                }
            }
        }

        void vector_factory_test() {
            auto &&comm = Comm::get_default();

            Vector vec(layout(comm, Traits::decide(), n_dofs), 0.2);
            {
                Read<Vector> read(vec);
                Range r = range(vec);

                for (int i = r.begin(); i != r.end(); ++i) {
                    utopia_test_assert(approxeq(vec.get(i), 0.2));
                }
            }
        }

        WrapperTest() = default;

    private:
        int n_dofs{100};
    };

    static void wrapper() {
#ifdef UTOPIA_WITH_PETSC
        WrapperTest<PetscMatrix, PetscVector>().run();
#endif

#ifdef WITH_BLAS
        WrapperTest<BlasMatrixd, BlasVectord>().run();
#endif  // WITH_BLAS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(wrapper);
}  // namespace utopia

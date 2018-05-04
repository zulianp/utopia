
#include "utopia_WrapperTest.hpp"
#include "utopia.hpp"

namespace utopia {

    template<class Matrix, class Vector, class Scalar>
    class WrapperTest {
    public:

        void run() {
            UTOPIA_RUN_TEST(vector_factory_test);
            UTOPIA_RUN_TEST(matrix_factory_test);
            UTOPIA_RUN_TEST(matrix_assembly_test);
        }

        void matrix_assembly_test() {
            Matrix mat = values(n_dofs, n_dofs, 0);

            //Assemble 1D laplacian with dirichlet nodes at the boundary
            {
                Write<Matrix> write(mat);
                Range rr = row_range(mat);

                //FIXME assumed row is owned by this proc completely.
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

            Vector expected = values(n_dofs, 0.);

            //Assemble expected result
            {
                Write<Vector> write(expected);
                Range r = range(expected);
                if (r.begin() == 0) {
                    expected.set(0, 10);
                }

                if (r.end() == expected.size().get(0)) {
                    expected.set(expected.size().get(0) - 1, 10);
                }
            }

            Vector vec = values(n_dofs, 10);
            assert(approxeq(expected, mat * vec));
        }

        void matrix_factory_test() {

            Matrix mat = values(n_dofs, n_dofs, 0.1);
            {
                Read<Matrix> read(mat);
                Range rr = row_range(mat);

                for (int i = rr.begin(); i != rr.end(); ++i) {
                    for (int j = 0; j != n_dofs; ++j) {
                        assert(approxeq(mat.get(i, j), 0.1));
                    }
                }
            }

            mat = identity(n_dofs, n_dofs);
            {
                Read<Matrix> read(mat);
                Range rr = row_range(mat);

                for (int i = rr.begin(); i != rr.end(); ++i) {
                    for (int j = 0; j != n_dofs; ++j) {
                        assert(approxeq(mat.get(i, j), i == j));
                    }
                }
            }
        }


        void vector_factory_test() {
            Vector vec = values(n_dofs, 0.2);
            {
                Read<Vector> read(vec);
                Range r = range(vec);

                for (int i = r.begin(); i != r.end(); ++i) {
                    assert(approxeq(vec.get(i), 0.2));
                }
            }
        }


        WrapperTest()
                : n_dofs(100) { }

    private:
        int n_dofs;

    };


    void runWrapperTest() {
        UTOPIA_UNIT_TEST_BEGIN("WrapperTest");
#ifdef WITH_PETSC
        WrapperTest<DMatrixd, DVectord, PetscScalar>().run();
#endif

#ifdef WITH_BLAS
        WrapperTest<Matrixd, Vectord, double>().run();
#endif //WITH_BLAS        
        UTOPIA_UNIT_TEST_END("WrapperTest");
    }
}


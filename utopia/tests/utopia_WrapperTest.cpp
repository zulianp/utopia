
#include "utopia_WrapperTest.hpp"
#include "utopia.hpp"

namespace utopia {

    template<class Matrix, class Vector, class Scalar>
    class WrapperTest {
    public:

        void run() {
            vectorTestFactory();
            matrixTestFactory();
            matrixAssemblyTest();
        }

        void matrixAssemblyTest() {
            Matrix mat = values(_n, _n, 0);

            //Assemble 1D laplacian with dirichlet nodes at the boundary
            {
                Write<Matrix> write(mat);
                Range rr = rowRange(mat);
                //Range cr = colRange(mat);

                //FIXME assumed row is owned by this proc completely.
                int rbegin = rr.begin();
                int rend = rr.end();

                if (rbegin == 0) {
                    mat.set(0, 0, 1);
                    rbegin += 1;
                }

                if (rend == mat.size().get(0)) {
                    mat.set(_n - 1, _n - 1, 1);
                    rend -= 1;
                }

                for (int i = rbegin; i != rend; ++i) {
                    mat.set(i, i, 2);
                    mat.set(i, i - 1, -1);
                    mat.set(i, i + 1, -1);
                }
            }

            Vector expected = values(_n, 0.);

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

            Vector vec = values(_n, 10);
            assert(approxeq(expected, mat * vec));
        }

        void matrixTestFactory() {

            Matrix mat = values(_n, _n, 0.1);
            {
                Read<Matrix> read(mat);
                Range rr = rowRange(mat);
                Range cr = colRange(mat);

                for (int i = rr.begin(); i != rr.end(); ++i) {
                    for (int j = cr.begin(); j != cr.end(); ++j) {
                        assert(approxeq(mat.get(i, j), 0.1));
                    }
                }
            }

            mat = identity(_n, _n);
            {
                Read<Matrix> read(mat);
                Range rr = rowRange(mat);
                Range cr = colRange(mat);

                for (int i = rr.begin(); i != rr.end(); ++i) {
                    for (int j = cr.begin(); j != cr.end(); ++j) {
                        assert(approxeq(mat.get(i, j), i == j));
                    }
                }
            }
        }


        void vectorTestFactory() {
            Vector vec = values(_n, 0.2);
            {
                Read<Vector> read(vec);
                Range r = range(vec);

                for (int i = r.begin(); i != r.end(); ++i) {
                    assert(approxeq(vec.get(i), 0.2));
                }
            }
        }


        WrapperTest()
                : _n(10) { }

    private:
        int _n;

    };


    void runWrapperTest() {
        std::cout << "Begin: WrapperTest" << std::endl;
#ifdef WITH_PETSC
        WrapperTest<DMatrixd, DVectord, PetscScalar> petsc;
        petsc.run();
#endif

#ifdef WITH_BLAS
        WrapperTest<Matrixd, Vectord, double> blas;
        blas.run();
        std::cout << "End:   WrapperTest" << std::endl;
#endif //WITH_BLAS        
    }
}


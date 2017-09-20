
#include "utopia_SpecTest.hpp"
#include "utopia.hpp"

#include <iostream>

namespace utopia {

    template<class Backend>
    class SpecTest {
    public:
        typedef typename Backend::Scalar Scalar;
        typedef typename Backend::Vector Vector;
        typedef typename Backend::Matrix Matrix;

        void run() {

            print_backend_info();
            /* Test vectors operation */
            UTOPIA_RUN_TEST(vec_ops_test);
            /* Test vector assignments */
            UTOPIA_RUN_TEST(vec_assign_test);
            /* Test matrix assignments */
            UTOPIA_RUN_TEST(matrix_assign_test);
            /* Test vector range */
            UTOPIA_RUN_TEST(vector_range_test);
            /* Test Matrix range */
            UTOPIA_RUN_TEST(matrix_range_test);



            
            /* Test matrix transpose assignment */
//            matrix_assigned_transposed_test();
//            /* Test assignment from range */
//            assign_from_range_test();
//            /* Test assignment to range */
//            testAssignToRange();
//            /* Test get size */
//            testSize();
//            /* BLAS 1 tests */
//            testScale();
//            testZaxpy();
//            testVectorReduce();
//            testDot();
//            testNorm2();
//            /* Test matrices operations */
//            testMatrixOps();
//            testMatrixScale();
//            testMatrixReduce();
//            testZaxpyMatrix();
//            testMatrixGEMV();
//            testMatrixGEMM();
        }

        void vec_ops_test() {
            testVecOp<Minus>(1, 2, -1);
            testVecOp<Plus>(1, 2, 3);
            testVecOp<EMultiplies>(2, 2, 4);
            // testVecOp<Divides>(2, 4, 0.5);
        }


        template<class Operation>
        void testVecOp(const Scalar lValue, const Scalar rValue, const Scalar result) {
            Backend &backend = Backend::Instance();
            Vector left, right, expected, actual;

            //BUILD FUNCTIONS
            backend.build(left, Size({_n, 1}), Values<Scalar>(lValue));
            backend.build(right, Size({_n, 1}), Values<Scalar>(rValue));
            backend.build(expected, Size({_n, 1}), Values<Scalar>(result));

            //OPERATION FUNCTIONS
            backend.apply(left, right, Operation(), actual);
            assert(backend.compare(expected, actual, ApproxEqual()));
        }

        void vec_assign_test() {
            const Scalar value = 1.;
            Backend &backend = Backend::Instance();
            Vector left_lvalue, right, right_rvalue;

            /* Build right vectors */
            backend.build(right, Size({_n, 1}), Values<Scalar>(value));
            backend.build(right_rvalue, Size({_n, 1}), Values<Scalar>(value));

            /* Execute reference assignment */
            backend.assign(left_lvalue, right);

            /* Test reference assignment */
            assert(backend.compare(left_lvalue, right, ApproxEqual()));


            /* Clear values in left_lvalue */
            backend.build(left_lvalue, Size({_n, 1}), Values<Scalar>(value + 1));

            /* Execute move assignment */
            backend.assign(left_lvalue, std::move(right_rvalue));

            /* Test result after assign with move */
            assert(backend.compare(left_lvalue, right, ApproxEqual()));
        }


        void matrix_assign_test() {
            const Scalar value = 1.;
            Backend &backend = Backend::Instance();
            Matrix left_lvalue, right, right_rvalue;

            /* Build right matrices */
            backend.build(right, Size({_n, _n}), Values<Scalar>(value));
            backend.build(right_rvalue, Size({_n, _n}), Values<Scalar>(value));

            /* Execute reference assignment */
            backend.assign(left_lvalue, right);

            /* Test reference assignment */
            assert(backend.compare(left_lvalue, right, ApproxEqual()));


            /* Clear values in left_lvalue */
            backend.build(left_lvalue, Size({_n, _n}), Values<Scalar>(value + 1));

            /* Execute move assignment */
            backend.assign(left_lvalue, std::move(right_rvalue));

            /* Test result after assign with move */
            assert(backend.compare(left_lvalue, right, ApproxEqual()));
        }

        void vector_range_test() {
            Backend &backend = Backend::Instance();
            Vector v0, v1;

            // Get number of MPI processes
            int comm_size = mpi_world_size();

            /* Build different vector sizes */
            backend.build(v0, Size({_n * comm_size, 1}), Zeros());
            backend.build(v1, Size({0, 1}), Zeros());

            /* Check ranges */
            if (backend.info().get_name() != "petsc") {
                assert(backend.range(v0).begin() == 0 &&
                       backend.range(v0).extent() == comm_size * _n &&
                       backend.range(v0).end() == comm_size * _n);
            } else {
                int rank = mpi_world_rank();
                assert(backend.range(v0).begin() == rank * _n &&
                       backend.range(v0).extent() == _n &&
                       backend.range(v0).end() == (rank + 1) * _n);
            }

            assert(backend.range(v1).begin() == 0 &&
                   backend.range(v1).extent() == 0 &&
                   backend.range(v1).end() == 0);
        }

        void matrix_range_test() {
            Backend &backend = Backend::Instance();
            Matrix m0, m1, m2;

            // Get number of MPI processes
            int comm_size = mpi_world_size();

            /* Build different matrices sizes */
            backend.build(m0, Size({_n * comm_size, _n * comm_size}), Zeros());
            backend.build(m1, Size({0, 0}), Zeros());
            backend.build(m2, Size({_n * comm_size, _n}), Zeros());

            /* Check ranges */
            if (backend.info().get_name() != "petsc") {
                assert(backend.rowRange(m0).begin() == 0 &&
                       backend.rowRange(m0).extent() == _n * comm_size &&
                       backend.rowRange(m0).end() == _n * comm_size);

                assert(backend.rowRange(m2).begin() == 0 &&
                       backend.rowRange(m2).extent() == _n * comm_size &&
                       backend.rowRange(m2).end() == _n * comm_size);

            } else {
                int rank = mpi_world_rank();

                assert(backend.rowRange(m0).begin() == rank * _n &&
                       backend.rowRange(m0).extent() == _n &&
                       backend.rowRange(m0).end() == (rank + 1) * _n);

                assert(backend.rowRange(m2).begin() == rank * _n &&
                       backend.rowRange(m2).extent() == _n &&
                       backend.rowRange(m2).end() == (rank + 1) * _n);
            }

            assert(backend.colRange(m0).begin() == 0 &&
                   backend.colRange(m0).extent() == _n * comm_size &&
                   backend.colRange(m0).end() == _n * comm_size);

            assert(backend.colRange(m2).begin() == 0 &&
                   backend.colRange(m2).extent() == _n &&
                   backend.colRange(m2).end() == _n);

            // Check empty matrix
            assert(backend.rowRange(m1).begin() == 0 &&
                   backend.rowRange(m1).extent() == 0 &&
                   backend.rowRange(m1).end() == 0);

            assert(backend.colRange(m1).begin() == 0 &&
                   backend.colRange(m1).extent() == 0 &&
                   backend.colRange(m1).end() == 0);

        }

        void matrix_assigned_transposed_test() {
            Backend &backend = Backend::Instance();
            Matrix left, right;

            /* Test transposed assignment for square matrix */
            backend.build(right, Size({_n, _n}), Zeros());
            backend.build(left, Size({_n, _n}), Zeros());

            /* Change entries of the matrix */
            {
                /* Acquiring lock for writing */
                backend.writeLock(right);

                Range rr = backend.rowRange(right);
                Range cr = backend.colRange(right);

                long rbegin = rr.begin();
                long rend = rr.end();

                long cbegin = cr.begin();
                long cend = cr.end();

                for (long j = rbegin; j < rend; ++j) {
                    for (long i = cbegin; i < cend; ++i) {
                        backend.set(right, j, i, i);
                    }
                }

                /* Releasing lock for writing */
                backend.writeUnlock(right);
            }

            /* Assign transposed matrix */
            backend.assignTransposed(left, right);

            /* Check result */
            {
                /* Acquiring lock to read */
                backend.readLock(left);

                Range rr = backend.rowRange(left);
                Range cr = backend.colRange(left);

                long rbegin = rr.begin();
                long rend = rr.end();

                long cbegin = cr.begin();
                long cend = cr.end();

                for (long j = rbegin; j < rend; ++j) {
                    for (long i = cbegin; i < cend; ++i) {
                        assert(backend.get(left, j, i) == j);
                    }
                }

                /* Releasing lock for reading */
                backend.readUnlock(left);
            }
        }

        void testAssignFromRange() {
            Backend &backend = Backend::Instance();

            /* Matrix test */
            Matrix left, right;

            /* Create test ranges */
            Range col_range(0, _n - 2);
            Range row_range(0, _n - 1);

            /* Test assignment using a range */
            backend.build(right, Size({_n - 1, _n - 2}), Values<Scalar>(1));
            backend.build(left, Size({_n, _n}), Zeros());


            backend.assignFromRange(left, right, row_range, col_range);

            /* Check operation result */
            {
                /* Acquiring lock to read */
                backend.readLock(left);

                Range rr = backend.rowRange(left);
                Range cr = backend.colRange(left);

                long rbegin = rr.begin();
                long rend = rr.end();

                long cbegin = cr.begin();
                long cend = cr.end();

                for (long j = rbegin; j < rend; ++j) {
                    for (long i = cbegin; i < cend; ++i) {
                        if (j < _n - 1 && i < _n - 2) {
                            assert(backend.get(left, j, i) == 1);
                        } else {
                            assert(backend.get(left, j, i) == 0);
                        }
                    }
                }

                /* Releasing lock for reading */
                backend.readUnlock(left);
            }

            /* Vector test */
            Vector left_v, right_v, expected;

            /* Vector assign range */
            Range v_range(1, _n - 1);
            Range c_range(0, 1);

            backend.build(right_v, Size({_n, 1}), Values<Scalar>(2));
            backend.build(expected, Size({_n - 2, 1}), Values<Scalar>(2));

            /* Assign vector from range */
            backend.assignFromRange(left_v, right_v, v_range, c_range);

            /* Check result */
            assert(backend.compare(left_v, expected, ApproxEqual()));

        }

        void testAssignToRange() {
            Backend &backend = Backend::Instance();

            /* Matrix test */
            Matrix left, right;

            /* Create test ranges */
            Range col_range(1, _n - 1);
            Range row_range(1, _n - 2);

            /* Test assignment using a range */
            backend.build(right, Size({_n - 2, _n - 1}), Values<Scalar>(1));
            backend.build(left, Size({_n, _n}), Zeros());

            if (backend.info().get_name() != "petsc") { //FIXME
                backend.assignToRange(left, right, row_range, col_range);

                /* Check operation result */
                {
                    /* Acquiring lock to read */
                    backend.readLock(left);

                    Range rr = backend.rowRange(left);
                    Range cr = backend.colRange(left);

                    long rbegin = rr.begin();
                    long rend = rr.end();

                    long cbegin = cr.begin();
                    long cend = cr.end();

                    for (long j = rbegin; j < rend; ++j) {
                        for (long i = cbegin; i < cend; ++i) {
                            if (j > 0 && j < _n - 2 && i > 0 && i < _n - 1) {
                                assert(backend.get(left, j, i) == 1);
                            } else {
                                assert(backend.get(left, j, i) == 0);
                            }
                        }
                    }

                    /* Releasing lock for reading */
                    backend.readUnlock(left);
                }
            }

            /* Vector test */
            Vector left_v, right_v;

            /* Vector assign range */
            Range v_range(1, _n - 1);
            Range c_range(0, 1);

            backend.build(left_v, Size({_n, 1}), Zeros());
            backend.build(right_v, Size({_n, 1}), Values<Scalar>(1));

            /* Assign vector from range */
            if (backend.info().get_name() != "petsc") { //FIXME
                backend.assignToRange(left_v, right_v, v_range, c_range);

                Range r = backend.range(left_v);

                /* Check result */
                for (   long j = r.begin();
                        j < r.end();
                        ++j) {
                    if (j > 0 && j < _n - 1) {
                        assert(backend.get(left_v, j) == 1);
                    } else {
                        assert(backend.get(left_v, j) == 0);
                    }
                }
            }
        }

        void testSize() {
            /* Test get size for matrices */
            Backend &backend = Backend::Instance();

            /* Matrix test */
            Matrix m1, m2, m3;
            /* Vector test */
            Vector v1, v2;

            /* Create different matrices */
            backend.build(m1, Size({_n, _n}), Zeros());
            backend.build(m2, Size({0, 0}), Zeros());
            backend.build(m3, Size({_n + 3, _n - 2}), Zeros());

            /* Create different vectors */
            backend.build(v1, Size({_n, 1}), Zeros());
            backend.build(v2, Size({0, 1}), Zeros());

            /* Test sizes */
            Size s;
            /* Test m1 */
            backend.size(m1, s);
            assert(s.nDims() == 2);
            assert(s.get(0) == _n && s.get(1) == _n);

            /* Test m2 */
            backend.size(m2, s);
            assert(s.nDims() == 2);
            assert(s.get(0) == 0 && s.get(1) == 0);

            /* Test m3 */
            backend.size(m3, s);
            assert(s.nDims() == 2);
            assert(s.get(0) == _n + 3 && s.get(1) == _n - 2);

            /* Test v1 */
            backend.size(v1, s);
            assert(s.nDims() == 1);
            assert(s.get(0) == _n);

            /* Test v2 */
            backend.size(v2, s);
            assert(s.nDims() == 1);
            assert(s.get(0) == 0);
        }

        /* BLAS 1 tests */

        void testScale() {
            /* Test get size for matrices */
            Backend &backend = Backend::Instance();

            /* Create vectors */
            Vector result, v1, v2, expected1, expected2;

            /* Build vectors */
            backend.build(v1, Size({_n, 1}), Values<Scalar>(2));
            backend.build(v2, Size({_n, 1}), Values<Scalar>(0));
            backend.build(expected1, Size({_n, 1}), Values<Scalar>(4));
            backend.build(expected2, Size({_n, 1}), Values<Scalar>(0));

            if (backend.info().get_name() != "petsc") { //FIXME
                /* Apply scaling */
                backend.scal(2, v1, result);

                /* Check result */
                assert(backend.compare(result, expected1, ApproxEqual()));

                /* Apply scaling on 0's vector */
                backend.scal(4, v2, result);

                /* Check result */
                assert(backend.compare(result, expected2, ApproxEqual()));
            }
        }

        void testZaxpy() {
            Backend &backend = Backend::Instance();

            /* Create vectors */
            Vector v1, v2, result, expected;

            backend.build(v1, Size({_n, 1}), Values<Scalar>(1));
            backend.build(v2, Size({_n, 1}), Values<Scalar>(3));
            backend.build(expected, Size({_n, 1}), Values<Scalar>(6));

            /* Compute result = s * v + result */
            backend.zaxpy(3, v1, v2, result);

            /* Check result */
            assert( backend.compare(result, expected, ApproxEqual()) );
        }

        void testVectorReduce() {
            Backend &backend = Backend::Instance();

            /* Create vectors */
            Vector v;

            backend.build(v, Size({_n, 1}), Values<Scalar>(2));

            /* Apply reduction */
            Scalar r = backend.reduce(v, Plus());

            assert(approxeq(r, 2 * _n));
        }

        void testMatrixReduce() {
            Backend &backend = Backend::Instance();

            /* Create matrix */
            Matrix m1, m2;
            backend.build(m1, Size({_n, _n}), Zeros());
            backend.build(m2, Size({_n, _n}), Values<Scalar>(2));

            if (backend.info().get_name() != "petsc") { //FIXME
                /* Apply reduction */
                Scalar m1_reduction = backend.reduce(m1, Plus());
                Scalar m2_reduction = backend.reduce(m2, Plus());

                /* Check result */
                assert(approxeq(m1_reduction, 0.0));
                assert(approxeq(m2_reduction, _n * _n * 2));
            }
        }

        void testDot() {
            Backend &backend = Backend::Instance();

            /* Create vectors */
            Vector v1, v2;

            backend.build(v1, Size({_n, 1}), Values<Scalar>(2));
            backend.build(v2, Size({_n, 1}), Values<Scalar>(3));

            Scalar d = backend.dot(v1, v2);

            assert(approxeq(d, 60.0));

        }

        void testNorm2() {
            Backend &backend = Backend::Instance();

            /* Create vectors */
            Vector v;

            backend.build(v, Size({4, 1}), Values<Scalar>(1));

            Scalar n2 = backend.norm2(v);

            assert(approxeq(n2, 2.0));
        }

        void testMatrixOps() {
            testMatrixOp<Plus>(1, 2, 3);
            testMatrixOp<Multiplies>(1, 1, _n);
            testMatrixOp<Minus>(1, 2, -1);
        }

        template <class Operation>
        void testMatrixOp(const Scalar lvalue, const Scalar rvalue, const Scalar result) {
            Backend &backend = Backend::Instance();

            Matrix left, right, res, expected;

            backend.build(left, Size({_n, _n}), Values<Scalar>(lvalue));
            backend.build(right, Size({_n, _n}), Values<Scalar>(rvalue));
            backend.build(expected, Size({_n, _n}), Values<Scalar>(result));

            backend.apply(left, right, Operation(), res);

            /* Check result */
            assert( backend.compare(res, expected, ApproxEqual()));
        }

        void testMatrixScale() {
            Backend &backend = Backend::Instance();

            Matrix m, res, expected;

            backend.build(m, Size({_n, _n}), Values<Scalar>(3));
            backend.build(expected, Size({_n, _n}), Values<Scalar>(9));

            if (backend.info().get_name() != "petsc") { //FIXME
                backend.scal(3, m, res);

                assert(backend.compare(res, expected, ApproxEqual()));
            }
        }

        void testMatrixGEMV() {
            Backend &backend = Backend::Instance();

            Matrix A;
            Vector x, y, expected;

            /* Create vectors and matrix */
            backend.build(A, Size({_n, _n}), Values<Scalar>(1));
            backend.build(x, Size({_n, 1}), Values<Scalar>(2));
            backend.build(y, Size({_n, 1}), Values<Scalar>(2));
            backend.build(expected, Size({_n, 1}), Values<Scalar>(4 * _n + 2));

            if (backend.info().get_name() != "petsc") { //FIXME
                /* Compute y = 2 * A * x + y */
                backend.gemv(2, A, x, 1, y);

                assert(backend.compare(y, expected, ApproxEqual()));
            }
        }

        void testMatrixGEMM() {
            Backend &backend = Backend::Instance();

            Matrix A, B, C, expected;

            if (backend.info().get_name() != "petsc") { //FIXME
                /* Create matrices */
                backend.build(A, Size({_n - 1, _n}), Values<Scalar>(1));
                backend.build(B, Size({_n, _n - 1}), Values<Scalar>(2));
                backend.build(C, Size({_n - 1, _n - 1}), Values<Scalar>(2));
                backend.build(expected, Size({_n - 1, _n - 1}), Values<Scalar>(4 * _n + 2));

                backend.gemm(2, A, B, 1, C);

                /* Check result */
                assert(backend.compare(C, expected, ApproxEqual()));

                /* Check other gemm versions */
                backend.build(C, Size({_n, _n}), Values<Scalar>(2));
                backend.build(expected, Size({_n, _n}), Values<Scalar>(4 * (_n - 1) + 2));
                backend.gemm(2, A, B, true, true, 1, C);

                /* Check result */
                assert(backend.compare(C, expected, ApproxEqual()));

                Vector v;
                backend.build(v, Size({_n - 1, 1}), Values<Scalar>(1));

                /* Check other gemm versions, result should be the same */
                backend.build(C, Size({1, _n - 1}), Values<Scalar>(2));
                backend.build(expected, Size({1, _n - 1}), Values<Scalar>(4 * (_n - 1) + 2));
                backend.gemm(2, v, B, true, true, 1, C);

                /* Check result */
                assert(backend.compare(C, expected, ApproxEqual()));
            }
        }

        void testZaxpyMatrix() {
            Backend &backend = Backend::Instance();

            Matrix left, right, result, expected;

            backend.build(left, Size({_n, _n}), Values<Scalar>(1));
            backend.build(right, Size({_n, _n}), Values<Scalar>(2));
            backend.build(expected, Size({_n, _n}), Values<Scalar>(4));

            backend.zaxpy(2, left, right, result);

            assert(backend.compare(result, expected, ApproxEqual()));

        }

        SpecTest()
                : _n(10) { }

        static void print_backend_info()
        {
            if(Utopia::Instance().verbose()) {
                std::cout << "\nBackend: " << Backend::Instance().info().get_name() << std::endl;
            }
        }      

    private:
        int _n;

    };


    void runSpecTest() 
    {
        UTOPIA_UNIT_TEST_BEGIN("SpecTest");
#ifdef WITH_BLAS        
        SpecTest<Backend<double, BLAS>>().run();
#endif //WITH_BLAS

#ifdef WITH_PETSC
        SpecTest<Backend<PetscScalar, PETSC>>().run();
#endif //WITH_PETSC
        UTOPIA_UNIT_TEST_END("SpecTest");
    }
}

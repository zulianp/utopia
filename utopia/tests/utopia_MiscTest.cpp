#include "utopia_Testing.hpp"

#include "utopia.hpp"
#include "utopia_CRSToBlockCRS.hpp"
#include "utopia_FactoryMethod.hpp"
#include "utopia_ILUDecompose.hpp"

namespace utopia {

    void crs_matrix() {
        std::vector<int> row_ptr = {0, 1, 3};
        std::vector<int> colidx = {0, 0, 1};
        std::vector<double> values = {1, -1, 1};

        CRSMatrix<std::vector<double>, std::vector<int>> mat(row_ptr, colidx, values, 2);
        CRSMatrix<std::vector<double>, std::vector<int>, 2> block_mat;
        convert(mat, block_mat);

        utopia_test_assert(3 == mat.nnz());
        utopia_test_assert(4 == block_mat.nnz());
    }

#ifdef UTOPIA_ENABLE_LAPACK
#ifdef UTOPIA_ENABLE_BLAS
    void test_lapack_eigen_solver() {
        using namespace utopia;
        using SizeType = Traits<BlasMatrixd>::SizeType;

        const int n = 9;
        BlasMatrixd A;
        A.dense(serial_layout(n, n), 0.0);
        std::vector<SizeType> is{0, 1, 2, 3, 4, 5, 6, 7, 8};

        {
            Write<BlasMatrixd> w_A(A);
            A.set_matrix(is,
                         is,
                         {+6.944444e+05,
                          0,
                          0,
                          0,
                          0,
                          0,
                          0,
                          0,
                          0,
                          0,
                          +1.388889e+06,
                          -5.555556e+04,
                          0,
                          -3.194444e+05,
                          -3.472222e+05,
                          0,
                          0,
                          0,
                          0,
                          -5.555556e+04,
                          +1.388890e+06,
                          0,
                          -3.472222e+05,
                          -3.194448e+05,
                          -5.555555e-02,
                          0,
                          -3.472222e-01,
                          0,
                          0,
                          0,
                          +6.944451e+05,
                          0,
                          0,
                          0,
                          0,
                          0,
                          0,
                          -3.194444e+05,
                          -3.472222e+05,
                          0,
                          +6.944444e+05,
                          -2.777778e+04,
                          0,
                          0,
                          0,
                          0,
                          -3.472222e+05,
                          -3.194448e+05,
                          0,
                          -2.777778e+04,
                          +6.944451e+05,
                          -3.472222e-01,
                          0,
                          -2.777778e-02,
                          0,
                          0,
                          -5.555555e-02,
                          0,
                          0,
                          -3.472222e-01,
                          +1.388889e+00,
                          0,
                          -3.194444e-01,
                          0,
                          0,
                          0,
                          0,
                          0,
                          0,
                          0,
                          +6.944444e-01,
                          0,
                          0,
                          0,
                          -3.472222e-01,
                          0,
                          0,
                          -2.777778e-02,
                          -3.194444e-01,
                          0,
                          +6.944444e-01});
        }

        BlasMatrixd B;
        B.dense(serial_layout(n, n), 0.0);
        {
            Write<BlasMatrixd> w_B(B);
            B.set_matrix(is, is, {+6.944444e+05, 0, 0, 0, 0, 0, 0, 0, 0, 0, +2.387253e+06, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  +2.387802e+06, 0, 0, 0, 0, 0, 0, 0, 0, 0, +6.944451e+05, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  +1.193627e+06, 0, 0, 0, 0, 0, 0, 0, 0, 0, +1.193901e+06, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  +1.841198e+00, 0, 0, 0, 0, 0, 0, 0, 0, 0, +6.944444e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  +9.205991e-01});
        }

        BlasMatrixd V;

        {
            Write<BlasMatrixd> w(A);
            A.set(0, 1, 2);
            A.set(1, 0, 2);
        }

        BlasVectord e;

        bool ok;
        ok = spd_geig_small(A, B, 0.3, e, V);
        utopia_test_assert(ok);
        ok = spd_geig(A, B, e, V);
        utopia_test_assert(ok);
        ok = spd_eig(A, e, V);
        utopia_test_assert(ok);
    }

    void factory_method() {
        FactoryMethod<BlasMatrixd> factory;
        auto ptr = factory.make();

        auto fact_1 = make_factory<BlasMatrixd, BlasMatrixd>();
        auto ptr_2 = fact_1->make();
    }

#endif  // UTOPIA_ENABLE_BLAS
#endif  // UTOPIA_ENABLE_LAPACK

    static void misc() {
        UTOPIA_RUN_TEST(crs_matrix);


#ifdef UTOPIA_ENABLE_LAPACK
#ifdef UTOPIA_ENABLE_BLAS
        UTOPIA_RUN_TEST(test_lapack_eigen_solver);
#endif
#endif
    }

    UTOPIA_REGISTER_TEST_FUNCTION(misc);
}  // namespace utopia

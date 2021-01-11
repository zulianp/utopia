#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_Vc.hpp"

#include "utopia_DeviceOperations.hpp"
#include "utopia_DeviceReduce.hpp"
#include "utopia_Views.hpp"

#include "utopia_vc_DeviceEigenDecomposition.hpp"

namespace utopia {

    class VcTest {
    public:
        using Scalar = double;
        // using Scalar = float;
        using vScalar = Vc::Vector<Scalar>;
        static const int N = vScalar::Size;

        Scalar ones[N] = {1, 1, 1, 1};
        Scalar twos[N] = {2, 2, 2, 2};

        Scalar scalar_test(int repeat, std::vector<Scalar> &in_data, std::vector<Scalar> &out_data) {
            Scalar ret = 0;
            const int data_size = in_data.size();

            for (int i = 0; i < repeat; ++i) {
                for (int k = 0; k < data_size; ++k) {
                    out_data[k] += std::sin(in_data[k]);
                }
            }

            for (int k = 0; k < data_size; ++k) {
                ret += out_data[k];
            }

            return ret;
        }

        Scalar vector_test(int repeat, std::vector<Scalar> &in_data, std::vector<Scalar> &out_data) {
            Scalar ret = 0;
            const int data_size = in_data.size();
            const int size_blocks = (data_size / vScalar::Size) * vScalar::Size;
            static const int simd_size = vScalar::Size;
            // const int santitized = data_size - size_blocks;

            for (int i = 0; i < repeat; ++i) {
                for (int k = 0; k < size_blocks; k += vScalar::Size) {
                    for (int d = 0; d < simd_size; ++d) {
                        out_data[k + d] += std::sin(in_data[k + d]);
                    }
                }

                for (int k = size_blocks; k < data_size; ++k) {
                    out_data[k] += std::sin(in_data[k]);
                }
            }

            for (int k = 0; k < data_size; ++k) {
                ret += out_data[k];
            }

            return ret;
        }

        // from unaligned to aligned
        Scalar simd_vector_test(int repeat, std::vector<Scalar> &in_data, std::vector<Scalar> &out_data) {
            Scalar ret = 0;
            const int data_size = in_data.size();
            const int size_blocks = (data_size / vScalar::Size) * vScalar::Size;

            for (int i = 0; i < repeat; ++i) {
                for (int k = 0; k < size_blocks; k += vScalar::Size) {
                    auto in = vScalar(&in_data[k], Vc::Unaligned);
                    auto out = vScalar(&out_data[k], Vc::Unaligned);

                    out += sin(in);
                    out.store(&out_data[k], Vc::Unaligned);
                }

                for (int k = size_blocks; k < data_size; ++k) {
                    out_data[k] += std::sin(in_data[k]);
                }
            }

            for (int k = 0; k < size_blocks; k += vScalar::Size) {
                auto out = vScalar(&out_data[k], Vc::Unaligned);
                ret += out.sum();
            }

            for (int k = size_blocks; k < data_size; ++k) {
                ret += out_data[k];
            }

            return ret;
        }

        Scalar simd_vector_aligned_test(int repeat, std::vector<Scalar> &u_in_data, std::vector<Scalar> &u_out_data) {
            using VectorA = std::vector<Scalar, Vc::Allocator<Scalar>>;

            Scalar ret = 0;
            const int data_size = u_in_data.size();
            const int size_blocks = (data_size / vScalar::Size) * vScalar::Size;

            VectorA in_data(data_size), out_data(data_size);

            for (int k = 0; k < data_size; ++k) {
                in_data[k] = u_in_data[k];
                out_data[k] = u_out_data[k];
            }

            for (int i = 0; i < repeat; ++i) {
                for (int k = 0; k < size_blocks; k += vScalar::Size) {
                    auto in = vScalar(&in_data[k], Vc::Aligned);
                    auto out = vScalar(&out_data[k], Vc::Aligned);

                    out += sin(in);
                    out.store(&out_data[k], Vc::Aligned);
                }

                for (int k = size_blocks; k < data_size; ++k) {
                    out_data[k] += std::sin(in_data[k]);
                }
            }

            for (int k = 0; k < size_blocks; k += vScalar::Size) {
                auto out = vScalar(&out_data[k], Vc::Aligned);
                ret += out.sum();
            }

            for (int k = size_blocks; k < data_size; ++k) {
                ret += out_data[k];
            }

            for (int k = 0; k < data_size; ++k) {
                u_out_data[k] = out_data[k];
            }

            return ret;
        }

        void compare() {
            int n = 300;

            std::cout << "SIMD size: " << N << std::endl;

            std::vector<Scalar> in_data(n, 1.0);
            std::vector<Scalar> out_data(n, 0);

            int repeat = 50000;

            Chrono c;
            c.start();
            Scalar ret1 = scalar_test(repeat, in_data, out_data);
            c.stop();
            double serial_time = c.get_seconds();

            std::cout << c << std::endl;

            //////////////////////////////////////////////////////////////////////////////

            std::fill(in_data.begin(), in_data.end(), 1.0);
            std::fill(out_data.begin(), out_data.end(), 0.0);

            c.start();
            Scalar ret2 = vector_test(repeat, in_data, out_data);
            c.stop();

            std::cout << c << std::endl;

            utopia_test_assert(approxeq(ret1, ret2, 1e-15));

            //////////////////////////////////////////////////////////////////////////////

            std::fill(in_data.begin(), in_data.end(), 1.0);
            std::fill(out_data.begin(), out_data.end(), 0.0);

            c.start();
            Scalar ret3 = simd_vector_test(repeat, in_data, out_data);
            c.stop();

            std::cout << c << std::endl;

            utopia_test_assert(approxeq(ret1, ret3, 1e-7));
            utopia_test_assert(c.get_seconds() < serial_time);

            //////////////////////////////////////////////////////////////////////////////

            std::fill(in_data.begin(), in_data.end(), 1.0);
            std::fill(out_data.begin(), out_data.end(), 0.0);

            c.start();
            Scalar ret4 = simd_vector_aligned_test(repeat, in_data, out_data);
            c.stop();

            std::cout << c << std::endl;

            utopia_test_assert(approxeq(ret1, ret4, 1e-7));
            utopia_test_assert(c.get_seconds() < serial_time);
        }

        void simd_with_vec() {
            StaticVector<vScalar, 2> vv1, vv2, vv3;

            vv1[0].load(ones, Vc::Unaligned);
            vv1[1].load(ones, Vc::Unaligned);

            vv2[0].load(twos, Vc::Unaligned);
            vv2[1].load(twos, Vc::Unaligned);

            vv3 = vv1 + vv2;

            Scalar result[N] = {-1, -1, -1, -1};

            for (int i = 0; i < 2; ++i) {
                vv3[i].store(result, Vc::Unaligned);
                for (int i = 0; i < N; ++i) {
                    utopia_test_assert(approxeq(3.0, result[i], 1e-7));
                }
            }
        }

        void simd_mat_vec() {
            StaticVector<vScalar, 2> x;
            StaticVector<vScalar, 2> y;
            StaticMatrix<vScalar, 2, 2> mat;

            x[0].load(ones, Vc::Unaligned);
            x[1].load(ones, Vc::Unaligned);

            mat(0, 0).load(twos, Vc::Unaligned);
            mat(0, 1).load(ones, Vc::Unaligned);
            mat(1, 0).load(ones, Vc::Unaligned);
            mat(1, 1).load(twos, Vc::Unaligned);

            y = mat * x;

            Scalar result[N] = {-1, -1, -1, -1};
            y[0].store(result, Vc::Unaligned);

            for (int i = 0; i < 2; ++i) {
                y[i].store(result, Vc::Unaligned);
                for (int i = 0; i < N; ++i) {
                    utopia_test_assert(approxeq(3.0, result[i], 1e-7));
                }
            }
        }

        void simd_mat_sum() {
            StaticMatrix<vScalar, 2, 2> mat;

            mat(0, 0).load(twos, Vc::Unaligned);
            mat(0, 1).load(ones, Vc::Unaligned);
            mat(1, 0).load(ones, Vc::Unaligned);
            mat(1, 1).load(twos, Vc::Unaligned);

            vScalar result = sum(mat);

            Scalar actual[N] = {-1, -1, -1, -1};
            Scalar expected[N] = {6, 6, 6, 6};
            result.store(actual, Vc::Unaligned);

            for (int i = 0; i < N; ++i) {
                utopia_test_assert(approxeq(expected[i], result[i], 1e-7));
            }
        }

        void simd_mat_det() {
            StaticMatrix<vScalar, 2, 2> mat;

            mat(0, 0).load(twos, Vc::Unaligned);
            mat(0, 1).load(ones, Vc::Unaligned);
            mat(1, 0).load(ones, Vc::Unaligned);
            mat(1, 1).load(twos, Vc::Unaligned);

            vScalar result = det(mat);

            Scalar actual[N] = {-1, -1, -1, -1};
            Scalar expected[N] = {3, 3, 3, 3};
            result.store(actual, Vc::Unaligned);

            for (int i = 0; i < N; ++i) {
                utopia_test_assert(approxeq(expected[i], result[i], 1e-7));
            }
        }

        template <class Expr, class Values>
        inline static void eig(const DeviceExpression<Expr> &expr,
                               Values &eigen_values,
                               StaticMatrix<vScalar, 2, 2> &eigen_vectors) {
            // if(eigen_values.size() == 2) {
            DeviceEigenDecompositionVc<Expr>::apply(expr.derived(), eigen_values, eigen_vectors);
            // } else {
            //     DeviceEigenDecomposition<Expr>::apply(expr.derived(), eigen_values, eigen_vectors);
            // }
        }

        void simd_mat_eigs() {
            StaticMatrix<vScalar, 2, 2> A, V;
            StaticVector<vScalar, 2> e;

            A(0, 0).load(twos, Vc::Unaligned);
            A(0, 1).load(ones, Vc::Unaligned);
            A(1, 0).load(ones, Vc::Unaligned);
            A(1, 1).load(twos, Vc::Unaligned);

            eig(A, e, V);

            for (int i = 0; i < 2; ++i) {
                std::cout << e[i] << std::endl;
            }

            // Scalar data[9] = {0.020000000000000,
            //                   0.000001250000000,
            //                   0.000001250000000,
            //                   0.000001250000000,
            //                   0.000004436491673,
            //                   0.0,
            //                   0.000001250000000,
            //                   0.0,
            //                   0.000004436491673};

            // StaticMatrix<vScalar, 3, 3> A, V;
            // StaticVector<vScalar, 3> e;

            // A(0, 0) = data[0];
            // A(0, 1) = data[1];
            // A(0, 2) = data[2];

            // A(1, 0) = data[3];
            // A(1, 1) = data[4];
            // A(1, 2) = data[5];

            // A(2, 0) = data[6];
            // A(2, 1) = data[7];
            // A(2, 2) = data[8];

            // vScalar result = det(A);

            // eig(A, e, V);

            // Scalar actual[N] = {-1, -1, -1, -1};
            // Scalar expected[N] = {3, 3, 3, 3};
            // result.store(actual, Vc::Unaligned);

            // for (int i = 0; i < N; ++i) {
            //     utopia_test_assert(approxeq(expected[i], result[i], 1e-7));
            // }
        }

        void run() {
            UTOPIA_RUN_TEST(simd_with_vec);
            UTOPIA_RUN_TEST(simd_mat_vec);
            UTOPIA_RUN_TEST(simd_mat_sum);
            UTOPIA_RUN_TEST(simd_mat_det);
            UTOPIA_RUN_TEST(simd_mat_eigs);
        }
    };

    static void Vc_ops() { VcTest().run(); }

    UTOPIA_REGISTER_TEST_FUNCTION(Vc_ops);
}  // namespace utopia
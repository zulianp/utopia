#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_Vc.hpp"

namespace utopia {

    class VcTest {
    public:
        // using Scalar = double;
        using Scalar = float;
        using vScalar = Vc::Vector<Scalar>;
        static const int N = vScalar::Size;

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
            // const int santitized = data_size - size_blocks;

            for (int i = 0; i < repeat; ++i) {
                for (int k = 0; k < size_blocks; k += vScalar::Size) {
                    for (int d = 0; d < vScalar::Size; ++d) {
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
            // const int santitized = data_size - size_blocks;

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

            // for (int k = 0; k < data_size; ++k) {
            //     ret += out_data[k];
            // }
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
            // const int santitized = data_size - size_blocks;

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

        void run() {
            // std::cout << vScalar::Size << std::endl;

            int n = 100;

            std::vector<Scalar> in_data(n, 1.0);
            std::vector<Scalar> out_data(n, 0);

            int repeat = 5000000;

            Chrono c;
            c.start();
            Scalar ret1 = scalar_test(repeat, in_data, out_data);
            c.stop();

            std::cout << c << std::endl;

            //////////////////////////////////////////////////////////////////////////////

            std::fill(in_data.begin(), in_data.end(), 1.0);
            std::fill(out_data.begin(), out_data.end(), 0.0);

            c.start();
            Scalar ret2 = vector_test(repeat, in_data, out_data);
            c.stop();

            std::cout << c << std::endl;

            // utopia_test_assert(approxeq(ret1, ret2, 1e-15));

            //////////////////////////////////////////////////////////////////////////////

            std::fill(in_data.begin(), in_data.end(), 1.0);
            std::fill(out_data.begin(), out_data.end(), 0.0);

            c.start();
            Scalar ret3 = simd_vector_test(repeat, in_data, out_data);
            c.stop();

            std::cout << c << std::endl;

            // utopia_test_assert(approxeq(ret1, ret3, 1e-7));

            //////////////////////////////////////////////////////////////////////////////

            std::fill(in_data.begin(), in_data.end(), 1.0);
            std::fill(out_data.begin(), out_data.end(), 0.0);

            c.start();
            Scalar ret4 = simd_vector_aligned_test(repeat, in_data, out_data);
            c.stop();

            std::cout << c << std::endl;

            // utopia_test_assert(approxeq(ret1, ret4, 1e-7));
        }
    };

    static void Vc_ops() { VcTest().run(); }

    UTOPIA_REGISTER_TEST_FUNCTION(Vc_ops);
}  // namespace utopia
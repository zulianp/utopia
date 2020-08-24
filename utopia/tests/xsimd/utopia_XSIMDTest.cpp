#include "utopia_Testing.hpp"
#include "utopia_xsimd.hpp"

namespace utopia {

    class XSIMDTest {
    public:
        template <int N>
        void test_vec() {
            double arr[N];
            std::fill(arr, arr + N, 1.0);

            host::Batch<double, N> b(arr);

            auto c = b + b;
            c.store_unaligned(arr);

            for (int i = 0; i < N; ++i) {
                std::cout << (arr[i]) << std::endl;
            }
        }

        void run() {
            test_vec<2>();
            // test_vec<3>();
        }
    };

    static void xsimd_ops() { XSIMDTest().run(); }

    UTOPIA_REGISTER_TEST_FUNCTION(xsimd_ops);
}  // namespace utopia

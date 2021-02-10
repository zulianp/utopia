#include "utopia_Base.hpp"

#include "utopia_simd_Assembler.hpp"
#include "utopia_simd_Assembler_v2.hpp"
#include "utopia_simd_Quadrature.hpp"

#include "utopia_ui.hpp"

// include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_Literal.hpp"

// #include "utopia_UniformHex8.hpp"

// std
#include <cmath>

using namespace utopia;

void simd_fe_v2(Input &in) {
    using Scalar = float;
    static constexpr int Dim = 3;
    static constexpr int Lanes = Vc::Vector<Scalar>::Size;
    using SIMDType = Vc::Vector<Scalar>;
    using Vector3 = simd_v2::Vector<Scalar, Dim>;
    // using Vector3 = simd_v1::Vector<SIMDType, Dim>;

    long repeat = 10000;
    in.get("repeat", repeat);
    Chrono c;

    /////////////////////////////////////////////////
    //// Standard Version
    /////////////////////////////////////////////////

    Scalar expected_value = 0.0;
    double elapsed_serial = 0.0;

    {
        Scalar v1[Lanes][Dim];
        Scalar v2[Lanes][Dim];
        Scalar dot_v[Lanes];

        for (int l = 0; l < Lanes; ++l) {
            dot_v[l] = 0.0;
            for (int i = 0; i < Dim; ++i) {
                v1[l][i] = Scalar(1.0);
                v2[l][i] = Scalar(3.0);
            }
        }

        c.start();
        for (long r = 0; r < repeat; ++r) {
            for (int l = 0; l < Lanes; ++l) {
                for (int i = 0; i < Dim; ++i) {
                    v1[l][i] *= Scalar(0.1);
                    v2[l][i] *= Scalar(0.1);
                    v1[l][i] += v2[l][i];
                }
            }

            for (int l = 0; l < Lanes; ++l) {
                for (int i = 0; i < Dim; ++i) {
                    dot_v[l] += v1[l][i] * v2[l][i];
                }
            }
        }

        for (int l = 0; l < Lanes; ++l) {
            expected_value += dot_v[l];
        }

        c.stop();
        elapsed_serial = c.get_seconds();
    }

    /////////////////////////////////////////////////
    //// SIMD Version
    /////////////////////////////////////////////////

    Scalar actual_value = 0.0;
    double elapsed_simd = 0;
    {
        Vector3 v1, v2;

        v1.set(Scalar(1.0));
        v2.set(Scalar(3.0));

        c.start();

        SIMDType dot_v = Scalar(0.);
        for (long r = 0; r < repeat; ++r) {
            v1 = v1 * Scalar(0.1);
            v2 = v2 * Scalar(0.1);
            v1 += v2;

            dot_v += dot(v1, v2);
        }

        actual_value = dot_v.sum();

        c.stop();
        elapsed_simd = c.get_seconds();
    }
    /////////////////////////////////////////////////

    // std::cout << expected_value << " == " << actual_value << std::endl;
    std::cout << elapsed_simd << " < " << elapsed_serial << " speed up " << (elapsed_serial / elapsed_simd)
              << std::endl;

    utopia_test_assert(approxeq(expected_value, actual_value, 1e-10));
    utopia_test_assert(elapsed_simd < elapsed_serial);

    // using SIMDType = Vc::Vector<double>;

    // int order = 1;
    // in.get("order", order);

    // simd_v1::Quadrature<SIMDType, 3> q_vec;

    // if (!simd_v1::Gauss<SIMDType>::Hex::get(order, q_vec)) {
    //     std::cerr << "[Error] could not find quadrature for order " << order << "\n";
    //     assert(false);
    //     return;
    // }

    // simd_v2::Vector<double, 3> h{1.0, 1.0, 1.0}, t{1.0, 1.0, 1.0};
    // simd_v2::Vector<SIMDType, 3> g, p_global;
    // UniformHex8<double> hex;
    // hex.set(t, h);

    // // RefHex8 hex;
    // // auto one = One<double>::value();
    // for (int qp = 0; qp < q_vec.n_points(); ++qp) {
    //     std::cout << "--------------\n";

    //     for (int i = 0; i < 8; ++i) {
    //         auto fi = hex.fun(i, q_vec.point(qp));
    //         std::cout << fi << std::endl;
    //     }

    //     std::cout << "--------------\n";
    // }

    // hex.grad(0, q_vec.point(0), g);
    // hex.point(q_vec.point(0), p_global);

    // disp(g);
    // disp(p_global);

    // PhysicalGradient<UniformHex8<double>, simd_v2::Quadrature<SIMDType, 3>> pg(q_vec);
    // Differential<UniformHex8<double>, simd_v2::Quadrature<SIMDType, 3>> dx(q_vec);

    // auto gp_e = pg.make(hex);
    // auto dx_e = dx.make(hex);

    // gp_e.get(0, 0, g);
    // disp(g);

    // std::cout << dx_e.get(0) << "\n";
}

UTOPIA_REGISTER_APP(simd_fe_v2);

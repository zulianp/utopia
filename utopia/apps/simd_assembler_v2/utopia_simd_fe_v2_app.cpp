#include "utopia_Base.hpp"

// #include "utopia_simd_Assembler.hpp"
#include "utopia_simd_Assembler_v2.hpp"
#include "utopia_simd_Quadrature_v2.hpp"

#include "utopia_ui.hpp"

// include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_Literal.hpp"

#include "utopia_simd_AssemblyView_v2.hpp"
#include "utopia_simd_UniformHex8.hpp"
#include "utopia_simd_UniformQuad4.hpp"

// std
#include <cmath>

using namespace utopia;

using Scalar = double;
// using Scalar = float;

static constexpr int Lanes = Vc::Vector<Scalar>::Size;
using Vector3 = simd_v2::Vector<Scalar, 3>;
using Vector2 = simd_v2::Vector<Scalar, 2>;

void simd_fe_v2(Input &in) {
    static constexpr int Dim = 3;
    using SIMDType = Vector3::SIMDType;

    long repeat = 10000;
    in.get("repeat", repeat);
    Chrono c;

    /////////////////////////////////////////////////
    //// Standard Version (This code should be also vectorizable by the compiler)
    /////////////////////////////////////////////////

    Scalar expected_value = 0.0;
    double elapsed_serial = 0.0;

    {
        // __attribute__((__aligned__(32))) Scalar v1[Dim][Lanes];
        // __attribute__((__aligned__(32))) Scalar v2[Dim][Lanes];
        // __attribute__((__aligned__(32))) Scalar dot_v[Lanes];

        Scalar v1[Dim][Lanes];
        Scalar v2[Dim][Lanes];
        Scalar dot_v[Lanes];

        for (int l = 0; l < Lanes; ++l) {
            dot_v[l] = 0.0;
        }

        for (int i = 0; i < Dim; ++i) {
            for (int l = 0; l < Lanes; ++l) {
                v1[i][l] = Scalar(1.0);
                v2[i][l] = Scalar(3.0);
            }
        }

        c.start();
        for (long r = 0; r < repeat; ++r) {
            for (int i = 0; i < Dim; ++i) {
                for (int l = 0; l < Lanes; ++l) {
                    v1[i][l] *= Scalar(0.1);
                    v2[i][l] *= Scalar(0.1);
                    v1[i][l] += v2[i][l];
                }
            }

            for (int i = 0; i < Dim; ++i) {
                for (int l = 0; l < Lanes; ++l) {
                    dot_v[l] += v1[i][l] * v2[i][l];
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

        SIMDType dot_v;
        dot_v = Scalar(0.);
        for (long r = 0; r < repeat; ++r) {
            v1 *= Scalar(0.1);
            v2 *= Scalar(0.1);
            v1 += v2;

            dot_v += dot(v1, v2);
        }

        actual_value = dot_v.sum();

        c.stop();
        elapsed_simd = c.get_seconds();
    }
    /////////////////////////////////////////////////

    std::cout << expected_value << " == " << actual_value << std::endl;
    std::cout << elapsed_simd << " < " << elapsed_serial << " speed up " << (elapsed_serial / elapsed_simd) << " x"
              << std::endl;

    utopia_test_assert(approxeq(expected_value, actual_value, 1e-10));
    utopia_test_assert(elapsed_simd < elapsed_serial);
}

UTOPIA_REGISTER_APP(simd_fe_v2);

void simd_fe_hex_v2(Input &in) {
    // static constexpr int Dim = 3;
    // using SIMDType = Vector3::SIMDType;

    int order = 1;
    in.get("order", order);

    simd_v2::Quadrature<Scalar, 3> q;

    if (!simd_v2::Gauss<Scalar>::Hex::get(order, q)) {
        std::cerr << "[Error] could not find quadrature for order " << order << "\n";
        assert(false);
        return;
    }

    Scalar h[3] = {1.0, 1.0, 1.0};
    Scalar t[3] = {1.0, 1.0, 1.0};

    simd_v2::UniformHex8<Scalar> hex;
    hex.set(t, h);

    for (int qp = 0; qp < q.n_points(); ++qp) {
        std::cout << "--------------\n";

        for (int i = 0; i < 8; ++i) {
            auto fi = hex.fun(i, q.point(qp));
            std::cout << fi << std::endl;
        }

        std::cout << "--------------\n";
    }

    auto p = q.point(0);
    // std::cout << "alignof: " << alignof(p) << " sizeof: " << sizeof(p) << std::endl;

    for (int i = 0; i < 3; ++i) {
        auto v = p[i] * h[i] + t[i];
        std::cout << (v) << std::endl;
    }

    Vector3 g, p_global;
    hex.grad(0, q.point(0), g);
    hex.point(q.point(0), p_global);

    std::cout << "G--------------\n";
    disp(g);
    std::cout << "P--------------\n";
    disp(p_global);
    std::cout << "--------------\n";

    // Using varying version here
    PhysicalGradient<simd_v2::UniformHex8<Scalar>, simd_v2::Quadrature<Scalar, 3>, Varying<>> pg(q);
    Differential<simd_v2::UniformHex8<Scalar>, simd_v2::Quadrature<Scalar, 3>> dx(q);

    auto gp_e = pg.make(hex);
    auto dx_e = dx.make(hex);

    gp_e.get(0, 0, g);
    disp(g);

    std::cout << dx_e.get(0) << "\n";
}

UTOPIA_REGISTER_APP(simd_fe_hex_v2);

void simd_fe_quad_v2(Input &in) {
    int order = 1;
    in.get("order", order);

    simd_v2::Quadrature<Scalar, 2> q;

    if (!simd_v2::Gauss<Scalar>::Quad::get(order, q)) {
        std::cerr << "[Error] could not find quadrature for order " << order << "\n";
        assert(false);
        return;
    }

    Scalar h[2] = {1.0, 1.0};
    Scalar t[2] = {1.0, 1.0};

    Vector2 g, p_global;
    simd_v2::UniformQuad4<Scalar> quad;
    quad.set(t, h);

    for (int qp = 0; qp < q.n_points(); ++qp) {
        std::cout << "--------------\n";

        for (int i = 0; i < 4; ++i) {
            auto fi = quad.fun(i, q.point(qp));
            std::cout << fi << std::endl;
        }

        std::cout << "--------------\n";
    }

    quad.grad(0, q.point(0), g);
    quad.point(q.point(0), p_global);

    disp(g);
    disp(p_global);

    // Using varying version here
    PhysicalGradient<simd_v2::UniformQuad4<Scalar>, simd_v2::Quadrature<Scalar, 2>, Varying<>> pg(q);
    Differential<simd_v2::UniformQuad4<Scalar>, simd_v2::Quadrature<Scalar, 2>> dx(q);

    auto gp_e = pg.make(quad);
    auto dx_e = dx.make(quad);

    gp_e.get(0, 0, g);
    disp(g);

    std::cout << dx_e.get(0) << "\n";
}

UTOPIA_REGISTER_APP(simd_fe_quad_v2);

void simd_fe_test(Input &in) {
    int order = 1;
    in.get("order", order);

    simd_v2::Quadrature<Scalar, 2> q;

    if (!simd_v2::Gauss<Scalar>::Quad::get(order, q)) {
        std::cerr << "[Error] could not find quadrature for order " << order << "\n";
        assert(false);
        return;
    }

    Scalar h[2] = {1.0, 1.0};
    Scalar t[2] = {1.0, 1.0};

    simd_v2::UniformQuad4<Scalar> quad;
    quad.set(t, h);

    // SIMDType p;

    for (int qp = 0; qp < q.n_points(); ++qp) {
        std::cout << "--------------\n";

        // std::cout << "HERE" << std::endl;
        simd_v2::Vector<Scalar, 2> &p = q.point(qp);
        // std::cout << "alignof: " << alignof(p) << " sizeof: " << sizeof(p) << std::endl;
        // std::cout << "THERE" << std::endl;

        // auto x = p[0];
        // std::cout << "alignof: " << alignof(x) << " sizeof: " << sizeof(x) << std::endl;

        for (int i = 0; i < 4; ++i) {
            auto fi = quad.fun(i, p);
            std::cout << fi << std::endl;
        }

        std::cout << "--------------\n";
    }
}

UTOPIA_REGISTER_APP(simd_fe_test);

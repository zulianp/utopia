#include "utopia_Base.hpp"
#include "utopia_ui.hpp"

// include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_Literal.hpp"
#include "utopia_simd_Quadrature.hpp"

// #include "utopia_UniformHex8.hpp"

#include "utopia_simd_Assembler_v2.hpp"

// std
#include <cmath>

using namespace utopia;

void simd_fe_v2(Input &in) {
    simd_v2::Vector<double, 3> v;
    v.set(1.0);

    v.scale(2.0);

    disp(v);

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


#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
//include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_Core.hpp"
#include "utopia_Jacobi.hpp"
#include "utopia_kokkos_Traits.hpp"
#include "utopia_trilinos_Poisson.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_StructuredGrid.hpp"
#include <cmath>

namespace utopia {

    static void sc_mesh()
    {
        using Mesh          = utopia::Grid2d;
        using SizeType      = Mesh::SizeType;
        using Elem          = Mesh::Elem;
        using Quadrature    = utopia::Quadrature<Elem, 2>;
        using FunctionSpace = utopia::FunctionSpace<Mesh, 1>;
        using Dev           = FunctionSpace::Device;

        using DevFunctionSpace = FunctionSpace::ViewDevice;
        using DofIndex         = DevFunctionSpace::DofIndex;


        TrilinosCommunicator world;

        if(world.size() > 2)
        {
            return;
        }

        std::array<SizeType, 2>  local_begin = { 0, 0 };
        std::array<SizeType, 2>  local_end = {1, 2};

        auto rank = world.rank();
        if(world.size() == 2) {

            local_begin[0] = 0;
            local_begin[1] = rank;

            local_end[0] = 1;
            local_end[1] = rank + 1;

        }

        Mesh mesh(world, { 1, 2 }, local_begin, local_end, {{1.0, 1.0}, {2.0, 2.0}});

        Quadrature quadrature;
        FunctionSpace space(mesh);
        auto space_view = space.view_device();
        auto q_view     = quadrature.view_device();

        Dev::parallel_for(space.local_element_range(), UTOPIA_LAMBDA(const SizeType &i) {
            Elem e;
            space_view.elem(i, e);

            DofIndex dofs;
            space_view.dofs(i, dofs);
        });


        // space.each_element(UTOPIA_LAMBDA(const SizeType &i, const Elem &elem) {
        //     // PhysicalPoint<Elem, Quadrature>    point(elem, q);
        //     PhysicalGradient<Elem, Quadrature> grad(elem, q);
        //     // ShapeFunction<Elem, Quadrature>    fun(elem, q);
        //     Differential<Elem, Quadrature>     dX(elem, q);
        //     ElementMatrix mat;

        //     //update(i, elem, grad, dX, ...)

        //     DofIndex dofs;

        //     mat.set(0.0);

        //     // Point p;
        //     // Grad g_trial, g_test;

        //     auto n = q.n_points();
        //     for(std::size_t k = 0; k < n; ++k) {
        //         for(std::size_t j = 0; j < elem.n_nodes(); ++j) {
        //             auto g_test = grad(j, k);
        //             mat(j, j) += dot(g_test, g_test) * dX(k);
        //             for(std::size_t l = j + 1; l < elem.n_nodes(); ++l) {
        //                 const auto v = dot(g_test, grad(l, k)) * dX(k);
        //                 mat(j, l) += v;
        //                 mat(l, j) += v;
        //             }
        //         }
        //     }

        //     space.dofs(i, dofs);
        //     //local 2 global here
        // });
    }

    UTOPIA_REGISTER_APP(sc_mesh);

    static void poisson_2D()
    {
        SizeType n = 20;
        Poisson<TpetraMatrix, TpetraVector> poisson(n);

        Chrono c;
        c.start();

        TpetraVector x = 0.0 * poisson.rhs();
        ConjugateGradient<TpetraMatrix, TpetraVector> cg;

        // auto prec = std::make_shared<PointJacobi<TpetraMatrix, TpetraVector>>();
        // prec->verbose(true);

        // auto prec = std::make_shared<Jacobi<TpetraMatrix, TpetraVector>>();
        // maybe put this in the preconditioner interface
        // prec->preconditioner_mode(true);
        // prec->max_it(50);
        auto prec = std::make_shared<InvDiagPreconditioner<TpetraMatrix, TpetraVector>>();

        cg.set_preconditioner(prec);
        cg.verbose(true);
        cg.max_it(n*n);
        cg.rtol(1e-8);
        cg.solve(poisson.laplacian(), poisson.rhs(), x);

        c.stop();

        std::cout << c << std::endl;

        // poisson.reinit();
        // write("x.m", x);
    }

    UTOPIA_REGISTER_APP(poisson_2D);
}

#endif //WITH_TRILINOS


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
        using SizeType      = Grid2d::SizeType;
        using Elem          = Grid2d::Elem;
        using QuadratureT   = Quadrature<Elem, 2>;
        using FunctionSpace = utopia::FunctionSpace<Grid2d, 1>;
        using DevFunctionSpace = FunctionSpace::DeviceView;
        using DofIndex         = DevFunctionSpace::DofIndex;
        using Dev = FunctionSpace::Device;

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

        Grid2d mesh(world, { 1, 2 }, local_begin, local_end, {{1.0, 1.0}, {2.0, 2.0}});

        // QuadratureT q;

        std::cout << mesh.local_element_range() << std::endl;

        FunctionSpace space(mesh);
        auto space_view = space.view_device();
        space.each_element(UTOPIA_LAMBDA(const SizeType &i, const Elem &elem) {
            DofIndex dofs;
            space_view.dofs(i, dofs);
            for(auto d : dofs) {
                std::cout << d << " ";
            }

            std::cout << "\n";
        });


        // Dev::parallel_for()


        // space.each_element(UTOPIA_LAMBDA(const SizeType &i, const Elem &elem) {
        //     // PhysicalPoint<Elem, QuadratureT>    point(elem, q);
        //     PhysicalGradient<Elem, QuadratureT> grad(elem, q);
        //     // ShapeFunction<Elem, QuadratureT>    fun(elem, q);
        //     Differential<Elem, QuadratureT>     dX(elem, q);
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

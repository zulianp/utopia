
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
#include "utopia_AssemblyView.hpp"
#include "utopia_Views.hpp"
#include "utopia_StencilView.hpp"
#include "utopia_GraphView.hpp"

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
        // using Point            = utopia::StaticVector<double, 2>;
        // using Grad             = utopia::StaticVector<double, 2>;
        using ElementMatrix    = utopia::StaticMatrix<double, 4, 4>;


        TrilinosCommunicator world;

        if(world.size() > 1)
        {
            return;
        }

        SizeType nx = 1;
        SizeType ny = 2;
        Mesh mesh(world, { nx, ny }, { 0, 0 }, { nx, ny }, {{ 1.0, 1.0 }, { 2.0, 2.0 }});

        Quadrature quadrature;
        FunctionSpace space(mesh);
        auto space_view = space.view_device();
        auto q_view     = quadrature.view_device();


        Graph<Mesh> graph;
        graph.init(mesh);

        // PhysicalPoint<FunctionSpace, Quadrature> points(space, quadrature);
        // auto p_view = points.view_device();

        //Host space processing
        PhysicalGradient<FunctionSpace, Quadrature> gradients(space, quadrature);

        //View extraction
        auto g_view = gradients.view_device();

        Differential<FunctionSpace, Quadrature> differentials(space, quadrature);
        auto d_view = differentials.view_device();

        Chrono c;
        c.start();

        SizeType n_assemblies = 0;
        SizeType * n_assemblies_ptr = &n_assemblies;

        double op_sum = 0;
        double * op_sum_ptr = &op_sum;

        Dev::parallel_for(
            space.local_element_range(),
            UTOPIA_LAMBDA(const SizeType &i)
        {
            Elem e;
            DofIndex dofs;
            ElementMatrix mat;

            space_view.elem(i, e);

            //element-wise extraction
            const auto grad = g_view.make(i, e);
            const auto dx   = d_view.make(i, e);

            mat.set(0.0);

            const auto n = grad.n_points();
            for(std::size_t k = 0; k < n; ++k) {
                for(std::size_t j = 0; j < grad.n_functions(); ++j) {
                    const auto g_test = grad(j, k);
                    mat(j, j) += dot(g_test, g_test) * dx(k);
                    for(std::size_t l = j + 1; l < grad.n_functions(); ++l) {
                        const auto v = dot(g_test, grad(l, k)) * dx(k);
                        mat(j, l) += v;
                        mat(l, j) += v;
                    }
                }
            }

            space_view.dofs(i, dofs);

            ++(*n_assemblies_ptr);
            (*op_sum_ptr) += norm1(mat);
        });

        c.stop();
        std::cout << n_assemblies << " assemblies " << c << std::endl;
        std::cout << op_sum << std::endl;
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

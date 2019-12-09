
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
        using IndexView1    = Grid2d::IndexView1;
        using ScalarView1   = Grid2d::ScalarView1;
        using SizeType      = Grid2d::SizeType;
        using NodeIndexView = Grid2d::NodeIndexView;
        using Elem          = Grid2d::Elem;
        using QuadratureT   = Quadrature<Elem, 2>;
        using FunctionSpace = utopia::FunctionSpace<Grid2d, 1>;
        using Point         = utopia::StaticVector<double, 2>;
        using Grad          = utopia::StaticVector<double, 2>;
        using ElementMatrix = utopia::StaticMatrix<double, 4, 4>;

        TrilinosCommunicator world;

        if(world.size() > 2)
        {
            return;
        }

        IndexView1  dims("dims", 2);
        IndexView1  local_begin("local_begin", 2);
        IndexView1  local_end("local_end", 2);
        Box<double, 2> box;

        box.min[0] = 1.0;
        box.min[1] = 1.0;

        box.max[0] = 2.0;
        box.max[1] = 2.0;

        dims[0] = 1;
        dims[1] = 2;

        auto rank = world.rank();
        if(world.size() == 2) {

            local_begin[0] = 0;
            local_begin[1] = rank;

            local_end[0] = 1;
            local_end[1] = rank + 1;

        } else {
            local_begin[0] = 0;
            local_begin[1] = 0;

            local_end[0] = 1;
            local_end[1] = 2;
        }

        Grid2d g(world, dims, local_begin, local_end, box);

        QuadratureT q;

        FunctionSpace space(g);
        space.each_element(UTOPIA_LAMBDA(const SizeType &i, const Elem &elem) {
            PhysicalPoint<Elem, QuadratureT>    point(elem, q);
            PhysicalGradient<Elem, QuadratureT> grad(elem, q);
            ShapeFunction<Elem, QuadratureT>    fun(elem, q);
            Differential<Elem, QuadratureT>     dX(elem, q);
            ElementMatrix mat;

            mat.set(0.0);

            Point p;
            Grad g_trial, g_test;

            auto n = point.size();
            for(std::size_t k = 0; k < n; ++k) {
                for(std::size_t j = 0; j < elem.n_nodes(); ++j) {
                    grad.get(j, k, g_test);

                    mat(j, j) += dot(g_test, g_test) * dX(k);

                    for(std::size_t l = j + 1; l < elem.n_nodes(); ++l) {
                        grad.get(l, k, g_trial);
                        auto v = dot(g_test, g_trial) * dX(k);
                        mat(j, l) += v;
                        mat(l, j) += v;
                    }
                }
            }
        });
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


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
        auto r = g.local_element_range();

        g.each_element(UTOPIA_LAMBDA(const SizeType &i, const Elem &elem) {
            ArrayView<double, Elem::NNodes> funs;

            elem.fun(Kokkos::subview(q.points(), 0, Kokkos::ALL()), funs);

            std::cout << world.rank() << ") " << i << std::endl;

            for(auto n : elem.nodes()) {
                std::cout << n << " ";
            }

            std::cout << "\n";
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

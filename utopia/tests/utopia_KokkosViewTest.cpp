
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
#include "utopia_Testing.hpp"

//include edsl components
#include "utopia_Core.hpp"
#include "utopia_Jacobi.hpp"
#include "utopia_kokkos_Traits.hpp"
#include "utopia_Views.hpp"
#include "test_problems/utopia_trilinos_Poisson3D.hpp"



#include <cmath>

namespace utopia {

    static void kokkos_vector_view()
    {
        using ViewType  = Kokkos::View<double *>;
        using ViewType2 = Kokkos::View<double **>;

        //utopia view wrapper forwards constructor arguments to view
        VectorView<ViewType> x("x", 10);
        x.set(1.0);
        x.axpy(2., x);

        x += 0.5 * x;
        x *= 2.0;

        const double x_dot_x = dot(x, x);
        // disp(x_dot_x);

        const SizeType n = 2;

        ViewType2 kokkos_x2("x2", n, 10);
        Kokkos::parallel_for(n, UTOPIA_LAMBDA(const int i ) {
            VectorView<ViewType> x2(Kokkos::subview(kokkos_x2, i, Kokkos::ALL()));
            x2.set(i);
        });

        // for(SizeType i = 0; i < n; ++i) {
        //     VectorView<ViewType> x2(Kokkos::subview(kokkos_x2, i, Kokkos::ALL()));
        //     x2.describe();
        // }
    }

    static void kokkos_matrix_view()
    {
        using ViewType  = Kokkos::View<double **>;
        using ViewType2 = Kokkos::View<double ***>;

        //kokkos view
        MatrixView<ViewType> A("A", 4, 4);

        A.set(1.0);
        A += 0.5 * A;

        SizeType n = 2;        
        ViewType2 kokkos_A2("A2", n, 4, 4);

        Kokkos::parallel_for(n, UTOPIA_LAMBDA(const int i ) {
            MatrixView<ViewType> A2(Kokkos::subview(kokkos_A2, i, Kokkos::ALL(), Kokkos::ALL()));
            A2.set(i);
            A2 += 0.5 * A2;
        });

        for(SizeType i = 0; i < n; ++i) {
            MatrixView<ViewType> A2(Kokkos::subview(kokkos_A2, i, Kokkos::ALL(), Kokkos::ALL()));
            A2.describe();
        }
    }

    static void kokkos_poisson_2D()
    {
        SizeType n = 20;

        Chrono c;
        c.start();
        Poisson<TpetraMatrix, TpetraVector> poisson(n);
        c.stop();

        std::cout << c << std::endl;

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

        write("x.m", x);
    }

    static void kokkos_view()
    {
        UTOPIA_RUN_TEST(kokkos_vector_view);
        UTOPIA_RUN_TEST(kokkos_matrix_view);
        UTOPIA_RUN_TEST(kokkos_poisson_2D);
    }

    UTOPIA_REGISTER_TEST_FUNCTION(kokkos_view);
}

#endif //WITH_TRILINOS

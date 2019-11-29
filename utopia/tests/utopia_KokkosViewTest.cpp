
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
#include "utopia_Testing.hpp"

//include edsl components
#include "utopia_Core.hpp"
#include "utopia_Jacobi.hpp"
#include "utopia_kokkos_Traits.hpp"
#include "utopia_Views.hpp"
#include "test_problems/utopia_trilinos_Poisson3D.hpp"
#include "utopia_trilinos.hpp"


#include <cmath>

#include <Tpetra_FECrsGraph.hpp>

namespace utopia {

    static void kokkos_vector_view()
    {
        using ViewType  = Kokkos::View<double *>;
        using ViewType2 = Kokkos::View<double **>;

        //utopia view wrapper forwards constructor arguments to view
        VectorView<ViewType> x(args__, "x", 10);
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
        MatrixView<ViewType> A(args__, "A", 4, 4);

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

    static void device_matrix_view()
    {
        using Dev      = utopia::Traits<TpetraVector>::Device;
        using SizeType = utopia::Traits<TpetraVector>::SizeType;

        SizeType n = 5;

        Poisson<TpetraMatrix, TpetraVector> poisson(n);

        TpetraMatrix L = poisson.laplacian();
        auto device_L = device_view(L);

        Dev::parallel_for(row_range(L), UTOPIA_LAMBDA(const SizeType &i) {
            device_L.atomic_add(i, i, 1.0);
        });
    }


    template<class Matrix>
    void assemble_periodic_laplacian_1D(Matrix &m)
    {
        // n x n matrix with maximum 3 entries x row
        Write<Matrix> w(m);
        Range r = row_range(m);
        auto n = size(m).get(0);

        for(SizeType i = r.begin(); i != r.end(); ++i) {
            if(i > 0) {
                m.set(i, i - 1, -1.0);
            } else {
                m.set(i, n-1, -1.0);
            }

            if(i < n-1) {
                m.set(i, i + 1, -1.0);
            } else {
                m.set(i, 0, -1.0);
            }

            if(i == 0 || i == n - 1) {
                m.set(i, i, 1.);
            } else {
                m.set(i, i, 2.0);
            }
        }

    }

    static void fe_crs_graph()
    {
        using Dev       = Traits<TpetraVector>::Device;
        using SizeType  = Traits<TpetraVector>::SizeType;
        using MapType   = Tpetra::Map<SizeType, SizeType>;
        using GraphType = Tpetra::FECrsGraph<SizeType, SizeType>;
        using View      = Kokkos::View<SizeType*>;

        Teuchos::RCP<const MapType> ownedRowMap, ownedPlusSharedRowMap;
        size_t maxNumEntriesPerRow = 3;

        TpetraMatrix mat;

        const SizeType rank = mat.comm().rank();
        const SizeType size = mat.comm().size();

        SizeType n_local = 10;
        SizeType n_global = mat.comm().sum(n_local);

        SizeType n_ghosts = (size>0) * 2;
        View index_list("il", n_local + n_ghosts);


        Dev::parallel_for(n_local, UTOPIA_LAMBDA(const SizeType &i) {
            index_list(i) = i + n_local * rank;
        });

        if(size > 0) {
            index_list(n_local)     = rank == 0? n_global -1 : n_local * rank - 1;
            index_list(n_local + 1) = rank == size-1? 0 : n_local * (rank + 1);
        }

        ownedRowMap = Teuchos::rcp(new MapType(n_global, n_local, 0, mat.comm().get()));
        ownedPlusSharedRowMap = Teuchos::rcp(
            new MapType(
                Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                index_list,
                0,
                mat.comm().get()
        ));

        auto graph =  Teuchos::rcp(new GraphType(
            ownedRowMap,
            ownedPlusSharedRowMap,
            maxNumEntriesPerRow
        ));

        graph->fillComplete();

        mat.raw_type().reset(
            new TpetraMatrix::crs_mat_type(graph)
        );

        // assemble_periodic_laplacian_1D(mat);
        // see example in /Users/zulianp/Desktop/code/installations/trilinos-git/packages/tpetra/core/example/Finite-Element-Assembly
        disp(mat);
    }

    static void kokkos_view()
    {
        UTOPIA_RUN_TEST(kokkos_vector_view);
        UTOPIA_RUN_TEST(kokkos_matrix_view);

        UTOPIA_RUN_TEST(device_matrix_view);
        UTOPIA_RUN_TEST(fe_crs_graph);

        if(mpi_world_size() == 1) {
            UTOPIA_RUN_TEST(kokkos_poisson_2D);
        }
    }

    UTOPIA_REGISTER_TEST_FUNCTION(kokkos_view);
}

#endif //WITH_TRILINOS

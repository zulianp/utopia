
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
#include "utopia_Testing.hpp"

#include "utopia_trilinos_FECrsGraph.hpp"
//include edsl components
#include "utopia_Core.hpp"
#include "utopia_Jacobi.hpp"
#include "utopia_kokkos_Traits.hpp"
#include "utopia_Views.hpp"
#include "utopia_trilinos.hpp"

#include <cmath>
#include <Tpetra_FECrsGraph.hpp>

namespace utopia {

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
            // A2.describe();
        }
    }

    static void device_matrix_view()
    {
        using Traits   = utopia::Traits<TpetraVector>;
        using Dev      = Traits::Device;
        using SizeType = Traits::SizeType;

        SizeType n = 5;

        TrilinosCommunicator comm;

        TpetraMatrix L; L.sparse(layout(comm, Traits::decide(), Traits::decide(), n, n), 3, 3);
        assemble_periodic_laplacian_1D(L);

        auto device_L = device_view(L);

        Dev::parallel_for(row_range(L), UTOPIA_LAMBDA(const SizeType &i) {
            device_L.atomic_add(i, i, 1.0);
        });
    }

    static void fe_crs_graph()
    {
        using Dev       = Traits<TpetraVector>::Device;
        using SizeType  = Traits<TpetraVector>::SizeType;
        using LocalSizeType = Traits<TpetraVector>::LocalSizeType;
        using MapType   = Tpetra::Map<LocalSizeType, SizeType>;
        using GraphType = Tpetra::FECrsGraph<LocalSizeType, SizeType>;
        using View      = Kokkos::View<SizeType*>;
        using DualView  = Kokkos::DualView<std::size_t*>;

        Teuchos::RCP<const MapType> ownedRowMap, ownedPlusSharedRowMap;
        // size_t maxNumEntriesPerRow = 3;

        TpetraMatrix mat;

        const SizeType rank = mat.comm().rank();
        const SizeType size = mat.comm().size();

        SizeType n_local = 10;
        SizeType n_global = mat.comm().sum(n_local);

        SizeType n_ghosts = (size>0) * 2;
        View index_list("il", n_local + n_ghosts);
        DualView nnz_z_row("nnz_z_row", n_local + n_ghosts);

        auto nnz_z_row_dev = nnz_z_row.view_device();

        Dev::parallel_for(n_local, UTOPIA_LAMBDA(const SizeType &i) {
            index_list(i) = i + n_local * rank;

            nnz_z_row_dev(i) = 3;

            if(i == 0 && size > 0) {
                index_list(n_local)     = rank == 0? n_global -1 : n_local * rank - 1;
                index_list(n_local + 1) = rank == size-1? 0 : n_local * (rank + 1);
                nnz_z_row_dev(n_local) = 3;
                nnz_z_row_dev(n_local + 1) = 3;
            }
        });

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
            // maxNumEntriesPerRow
            nnz_z_row
        ));

        Range r(ownedRowMap->getMinGlobalIndex(), ownedRowMap->getMaxGlobalIndex() + 1);
        SizeType cols[3];
        for(SizeType i = r.begin(); i != r.end(); ++i) {
            cols[0] = i;

            if(i > 0) {
                cols[1] = i-1;
            } else {
                cols[1] = n_global-1;
            }

            if(i < n_global-1) {
                cols[2] = i + 1;
            } else {
                cols[2] = 0;
            }

            std::sort(std::begin(cols), std::end(cols));
            graph->insertGlobalIndices(i, 3, cols);
        }

        graph->fillComplete();

        mat.raw_type().reset(
            new TpetraMatrix::CrsMatrixType(graph)
        );

        assemble_periodic_laplacian_1D(mat);
        // see example in /Users/zulianp/Desktop/code/installations/trilinos-git/packages/tpetra/core/example/Finite-Element-Assembly
        // disp(mat);

        assemble_periodic_laplacian_1D(mat);
    }

    static void kokkos_view()
    {
        UTOPIA_RUN_TEST(kokkos_vector_view);
        UTOPIA_RUN_TEST(kokkos_matrix_view);

        UTOPIA_RUN_TEST(device_matrix_view);
        UTOPIA_RUN_TEST(fe_crs_graph);
    }

    UTOPIA_REGISTER_TEST_FUNCTION(kokkos_view);
}

#endif //WITH_TRILINOS

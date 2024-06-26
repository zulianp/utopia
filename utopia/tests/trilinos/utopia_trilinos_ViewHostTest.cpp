
#include "utopia_Base.hpp"

#ifdef UTOPIA_ENABLE_TRILINOS
#include "utopia_Testing.hpp"

#include "utopia_trilinos_FECrsGraph.hpp"
// include edsl components
#include "utopia_Core.hpp"
#include "utopia_Device.hpp"
#include "utopia_Jacobi.hpp"
#include "utopia_Views.hpp"

#include "utopia_kokkos_Traits.hpp"
#include "utopia_trilinos.hpp"

#include "utopia_trilinos_ViewHost.hpp"

namespace utopia {

    static void vector_view_host() {
        auto &&comm = TrilinosCommunicator::get_default();
        TpetraVector vec(layout(comm, 1, comm.size()), 1);

        {
            // LocalViewHost<TpetraVector> v(vec);
            auto v = local_view_host(vec);
            v.set(0, comm.rank());
        }

        {
            // LocalViewHost<TpetraVector> v(vec);
            auto v = local_view_host(vec);
            for (int r = 0; r < comm.size(); r++) {
                if (r == comm.rank()) {
                    if (r != v.get(0)) {
                        utopia::out() << v.get(0) << "\n";
                        utopia_test_assert(false);
                    }
                }

                comm.barrier();
            }
        }
    }

    // static void matrix_view_host() {
    //     auto &&comm = TrilinosCommunicator::get_default();
    //     TpetraMatrix mat;

    //     auto vl = layout(comm, 1, comm.size());
    //     mat.sparse(square_matrix_layout(vl), 2, 2);

    //     {
    //         Write<TpetraMatrixd> w(mat);

    //         auto rr = row_range(mat);
    //         for (auto i = rr.begin(); i < rr.end(); i++) {
    //             mat.set(i, i, i);
    //         }
    //     }

    //     ViewHost<TpetraMatrix> v(mat);

    //     auto rr = row_range(mat);
    //     for (auto i = rr.begin(); i < rr.end(); i++) {
    //         v.atomic_add(i, i, 1);
    //     }
    // }

    static void kokkos_view_host() {
        UTOPIA_RUN_TEST(vector_view_host);
        // UTOPIA_RUN_TEST(matrix_view_host);
    }

    UTOPIA_REGISTER_TEST_FUNCTION(kokkos_view_host);
}  // namespace utopia

#endif  // UTOPIA_ENABLE_TRILINOS

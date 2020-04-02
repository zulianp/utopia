#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_assemble_laplacian_1D.hpp"
#include "utopia_ParallelTestRunner.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class SubCommAlgebraTest final : public AlgebraUnitTest<Vector> {
    public:
        using Traits       = utopia::Traits<Vector>;
        using Scalar       = typename Traits::Scalar;
        using SizeType     = typename Traits::SizeType;
        using IndexSet     = typename Traits::IndexSet;
        using Comm         = typename Traits::Communicator;
        using Layout       = typename Traits::Layout;
        using MatrixLayout = typename Traits::MatrixLayout;

        void sum_vectors()
        {
            Vector v1(layout(this->comm(), 2, this->comm().size() * 2), 1.0);
            Vector v2(layout(v1), 2.0);
            Vector v3 = v1 + v2;
            Scalar norm_v3 = norm1(v3);
            utopia_test_assert(approxeq(norm_v3, this->comm().size() * 3 * 2, device::epsilon<Scalar>()));
        }

        void run()
        {
            sum_vectors();
        }

    };

    void sub_comm_algebra()
    {
#ifdef WITH_PETSC
        run_parallel_test< SubCommAlgebraTest<PetscMatrix, PetscVector> >();
#endif //WITH_PETSC

#ifdef WITH_TRILINOS
        run_parallel_test< SubCommAlgebraTest<TpetraMatrix, TpetraVector> >();
#endif //WITH_PETSC
    }

    UTOPIA_REGISTER_TEST_FUNCTION(sub_comm_algebra);

}
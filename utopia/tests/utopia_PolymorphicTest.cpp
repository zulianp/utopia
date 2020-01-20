#include "utopia_Testing.hpp"
#include "utopia.hpp"
#include <string>
#include <cassert>

#include "utopia_DynamicTypeDistributedMatrix.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class PolymorphicTest final {
    public:
        using Scalar   = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;

        //base classes
        using DistributedMatrix = utopia::DistributedMatrix<Scalar, SizeType>;
        using DistributedVector = utopia::DistributedVector<Scalar, SizeType>;

        void convenience_wrapper()
        {
            std::shared_ptr<DistributedMatrix> mat = std::make_shared<Matrix>();
            std::shared_ptr<DistributedVector> vec = std::make_shared<Vector>();

            DynamicTypeDistributedMatrix<Scalar, SizeType, Matrix, Vector> * mat_dyn;
        }

        void run()
        {
           UTOPIA_RUN_TEST(convenience_wrapper);
        }
    };

    static void polymorphic()
    {
#ifdef WITH_PETSC
    PolymorphicTest<PetscMatrix, PetscVector>().run();
#endif //WITH_PETSC
    }

    UTOPIA_REGISTER_TEST_FUNCTION(polymorphic);
}

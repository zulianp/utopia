#include "utopia_Testing.hpp"
#include "utopia.hpp"
#include <string>
#include <cassert>

#include "utopia_petsc.hpp"
#include "utopia_AbstractVector.hpp"

namespace utopia {

#ifdef WITH_PETSC

    class PolymorphicTest final {
    public:
        using Scalar   = typename Traits<PetscVector>::Scalar;
        using SizeType = typename Traits<PetscVector>::SizeType;

        //base classes
        // using DistributedMatrix = utopia::DistributedMatrix<Scalar, SizeType>;
        using AbstractVector = utopia::AbstractVector<Scalar, SizeType>;

        void convenience_wrapper()
        {
            const SizeType n = 10;

            // std::shared_ptr<DistributedMatrix> mat    = std::make_shared<Matrix>(local_identity(n, n));
            // std::shared_ptr<DistributedVector> vec    = std::make_shared<Vector>(local_values(n, 2.0));
            // std::shared_ptr<DistributedVector> result = std::make_shared<Vector>(local_zeros(n));
            // mat->multiply(*vec, *result);

            std::shared_ptr<AbstractVector> v = std::make_shared<Wrapper<PetscVector>>( local_values(n, 2.0) );
#ifdef WITH_TRILINOS
#ifdef UTOPIA_TPETRA_SIZE_TYPE
            v = std::make_shared<Wrapper<TpetraVector>>( local_values(n, 2.0) );
#endif //UTOPIA_TPETRA_SIZE_TYPE
#endif //WITH_TRILINOS

        }

        void run()
        {
           UTOPIA_RUN_TEST(convenience_wrapper);
        }
    };

    static void polymorphic()
    {
        PolymorphicTest().run();

    }

    UTOPIA_REGISTER_TEST_FUNCTION(polymorphic);

#endif //WITH_PETSC
}

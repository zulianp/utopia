#include "utopia_Testing.hpp"
#include "utopia.hpp"
#include <string>
#include <cassert>

#include "utopia_petsc.hpp"
#include "utopia_AbstractVector.hpp"
#include "utopia_ObjectFactory.hpp"

namespace utopia {

//Register types
#ifdef WITH_PETSC

    UTOPIA_FACTORY_REGISTER_VECTOR(PetscVector);

#ifdef WITH_TRILINOS
    UTOPIA_FACTORY_REGISTER_VECTOR(TpetraVector);
#endif

    class PolymorphicTest final {
    public:
        using Scalar   = typename Traits<PetscVector>::Scalar;
        using SizeType = typename Traits<PetscVector>::SizeType;

        using DefaultFactory = utopia::AlgebraFactory<Scalar, SizeType>;

        //base classes
        // using DistributedMatrix = utopia::DistributedMatrix<Scalar, SizeType>;
        using AbstractVector = utopia::AbstractVector<Scalar, SizeType>;

        void convenience_wrapper()
        {
            const SizeType n = 10;

#ifdef WITH_TRILINOS
#ifdef UTOPIA_TPETRA_SIZE_TYPE
            //if types are the same TirlinosFactory == DefaultFactory
            InputParameters params;
            params.set("default-backend", "trilinos");
            DefaultFactory::init(params);

#endif //UTOPIA_TPETRA_SIZE_TYPE
#endif //WITH_TRILINOS


            auto v = DefaultFactory::new_vector();
            v->values(n, 2.0);
            v->describe();
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

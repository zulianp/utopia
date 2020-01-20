#include "utopia_Testing.hpp"
#include "utopia.hpp"
#include <string>
#include <cassert>

namespace utopia {

    template<
        typename Scalar,
        typename SizeType,
        class ConcreteMatrixType,
        class ConcreteVectorType
    >
    class DynamicTypeDistributedMatrix :
        public DistributedMatrix<Scalar, SizeType>,
        public BLAS2Matrix<DistributedMatrix<Scalar, SizeType>,
                           DistributedVector<Scalar, SizeType>> {
    public:

        virtual ~DynamicTypeDistributedMatrix() {}

        //route abstract type methods to concrete methods


    };


    template<class Matrix, class Vector>
    class PolymorphicTest final {
    public:
        using Scalar   = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;

        void run()
        {
            std::shared_ptr< DistributedMatrix<Scalar, SizeType> > mat = std::make_shared<Matrix>();
            std::shared_ptr< DistributedVector<Scalar, SizeType> > vec = std::make_shared<Vector>();
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

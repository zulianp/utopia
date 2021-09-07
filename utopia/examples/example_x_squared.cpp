#include "utopia.hpp"

namespace utopia {
    template <class Vector>
    class XSquared : public FunctionBase<Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        XSquared();
        ~XSquared();
        bool value(const Vector &x, Scalar);
        bool gradient(const Vector &x, const Vector &g);
    };
}  // namespace utopia

int main(int argc, char **argv) {
    using namespace utopia;

#ifdef WITH_PETSC
    using MatrixT = PetscMatrix;
    using VectorT = PetscVector;
#else
#ifdef WITH_TRILINOS
    using MatrixT = TpetraMatrixd;
    using VectorT = TpetraVectord;
#endif
#endif

    Utopia::Init(argc, argv);
    std::cout << "Working \n";
    return Utopia::Finalize();
}  // namespace )
#include "utopia.hpp"
#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Version.hpp"

namespace utopia {
    template <class Vector>
    class XSquared final : public FunctionBase<Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        XSquared() = default;

        bool value(const Vector &x, Scalar value) {
            if (x.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            value = sum(x);

            return true;
        }
        bool gradient(const Vector &x, const Vector &g) {
            if (x.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }
            g = 2 * x;

            return true;
        }
    };

    class NewGradientDescent {
    public:
        NewGradientDescent() : lr_(0.1), tol_(0.1) { std::cout << "alright\n"; }

    private:
        double lr_;
        double tol_;
        // XSquared<class Vector> func;
    };

}  // namespace utopia

template <class Vector>
void test() {
    using namespace utopia;
    std::cout << "Running test..." << std::endl;

    using Comm = typename Traits<Vector>::Communicator;
    using SizeType = typename Traits<Vector>::SizeType;
    using Scalar = typename Traits<Vector>::Scalar;

    Comm comm = Comm::get_default();
    SizeType n_local = 10;
    SizeType n_global = n_local * comm.size();
    auto l = layout(comm, n_local, n_global);

    Vector x, g;

    x.values(l, 10);
    g.values(l, 0);

    XSquared<Vector> newXSquared;
}

int main(int argc, char **argv) {
    using namespace utopia;

// Checking which backend is active.
#ifdef UTOPIA_WITH_PETSC
    using MatrixT = PetscMatrix;
    using VectorT = PetscVector;
#endif

    Utopia::Init(argc, argv);

    test<VectorT>();

    return Utopia::Finalize();
}  // namespace )

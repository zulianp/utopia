#include <utopia.hpp>

template <class Matrix, class Vector>
class Rosenbrock2DFunction : public utopia::Function<Matrix, Vector> {
public:
    typedef UTOPIA_SCALAR(Matrix) Scalar;

    Rosenbrock2DFunction() {}

    bool value(const Vector &point, Scalar &result) const override {
        using namespace utopia;

        assert(point.size() == 2);
        assert(!utopia::is_parallel<Matrix>::value ||
               point.comm().size() == 1 && "does not work for parallel matrices");

        const Read<Vector> read(point);

        const Scalar x = point.get(0);
        const Scalar y = point.get(1);

        result = 100.0 * pow(y - x * x, 2.0) + pow(1.0 - x, 2.0);
        return true;
    }

    bool gradient(const Vector &point, Vector &result) const override {
        using namespace utopia;

        assert(point.size() == 2);
        assert(!utopia::is_parallel<Matrix>::value ||
               point.comm().size() == 1 && "does not work for parallel matrices");

        if (empty(result)) {
            result.zeros(layout(point));
        }

        const Read<Vector> read(point);
        const Write<Vector> write(result);

        const Scalar x = point.get(0);
        const Scalar y = point.get(1);

        result.set(0, (400.0 * x * x * x - 400 * x * y + 2.0 * x - 2.0));
        result.set(1, 200.0 * (y - x * x));
        return true;
    }

    bool hessian(const Vector &point, Matrix &result) const override {
        using namespace utopia;

        assert(point.size() == 2);
        assert(!utopia::is_parallel<Matrix>::value ||
               point.comm().size() == 1 && "does not work for parallel matrices");

        if (empty(result)) {
            result.dense(square_matrix_layout(layout(point)));
        }

        const Read<Vector> read(point);
        const Write<Matrix> write(result);

        const Scalar x = point.get(0);
        const Scalar y = point.get(1);
        const Scalar mixed = -400.0 * x;

        result.set(0, 0, 1200 * x * x - 400 * y + 2);
        result.set(0, 1, mixed);
        result.set(1, 0, mixed);
        result.set(1, 1, 200.0);
        return true;
    }
};

int main(int argc, char **argv) {
    using namespace utopia;

#ifdef WITH_PETSC
    using MatrixT = PetscMatrix;
    using VectorT = PetscVector;
#else
    using MatrixT = BlasMatrixd;
    using VectorT = BlasVectord;
#endif
    Utopia::Init(argc, argv);

    InputParameters params;
    params.init(argc, argv);

    {  // Only in the main file: put this scope so that the petsc objects will be destroyed before the call to finalize

        // instatiating Rosenbrock 2D banana function
        Rosenbrock2DFunction<MatrixT, VectorT> rosenbrock_fun;

        VectorT rosenbrock_exact, x;

        // exact solution to our problem
        rosenbrock_exact.values(serial_layout(2), 1);

        // constructing initial guess
        x.values(serial_layout(2), 2);

        // nonlinear solve
        Newton<MatrixT, VectorT> newton;
        newton.solve(rosenbrock_fun, x);

        // comparing obtained solution with exact one

        std::stringstream ss;
        ss << "Correct solution: " << (approxeq(x, rosenbrock_exact) ? "true." : "false.") << std::endl;

        x.comm().root_print(ss.str());
    }

    return Utopia::Finalize();
}

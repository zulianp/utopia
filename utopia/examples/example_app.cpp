#include "utopia_Main.hpp"

#ifdef UTOPIA_WITH_PETSC
using Matrix_t = utopia::PetscMatrix;
using Vector_t = utopia::PetscVector;
#else
#ifdef UTOPIA_WITH_TRILINOS
using Matrix_t = utopia::TpetraMatrixd;
using Vector_t = utopia::TpetraVectord;
#else
using Matrix_t = utopia::BlasMatrixd;
using Vector_t = utopia::BlasVectord;
#endif
#endif

// ./examples/example_app -app solve -A <path_to_A> -b <path_to_b> &&  ./examples/example_app -app print
int main(const int argc, char *argv[]) { return UTOPIA_MAIN(argc, argv); }

void solve(utopia::Input &in) {
    using namespace utopia;
    using LinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;

    Path path_A = "A.bin", path_b = "b.bin", path_x = "out.bin";
    in.get("A", path_A);
    in.get("b", path_b);
    in.get("x", path_x);

    Matrix_t A;
    Vector_t b, x;

    if (!read(path_A, A) || !read(path_b, b)) {
        err() << "Unable to read input\n";
        return;
    }

    assert(!b.empty());
    assert(!A.empty());

    x.zeros(layout(b));

    LinearSolver_t solver;
    solver.read(in);

    solver.solve(A, b, x);

    if (!write(path_x, x)) {
        err() << "Unable to write output\n";
    }
}

UTOPIA_REGISTER_APP(solve);

void print(utopia::Input &in) {
    using namespace utopia;

    Path path_x = "out.bin";
    in.get("x", path_x);

    Vector_t x;
    if (read(path_x, x)) {
        disp(x);
    } else {
        err() << "Unable to read input\n";
    }
}

UTOPIA_REGISTER_APP(print);

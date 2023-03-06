#include "utopia_Main.hpp"

#ifdef UTOPIA_ENABLE_PETSC
using Matrix_t = utopia::PetscMatrix;
using Vector_t = utopia::PetscVector;
#else
#ifdef UTOPIA_ENABLE_TRILINOS
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

    Matrix_t A;
    Vector_t b, x;
    LinearSolver_t solver;
    Path path_A = "A.bin", path_b = "b.bin", path_x = "out.bin";

    Options opts;  // ./example_app -app solve  --show
    opts.add_option("A", path_A, "Path to system matrix (.bin or .mm)")
        .add_option("b", path_b, "Path to the right-hand side")
        .add_option("x", path_x, "Path to the output file")
        .add_option("linear_solver",
                    solver,
                    "Customize the solver (only works with JSON config file)\ne.g. @file settings.json");

    if (!opts.parse(in)) {
        return;
    }

    if (!read(path_A, A) || !read(path_b, b)) {
        err() << "Unable to read input\n";
        return;
    }

    assert(!b.empty());
    assert(!A.empty());

    x.zeros(layout(b));
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

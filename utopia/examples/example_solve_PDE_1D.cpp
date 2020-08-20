#include "utopia.hpp"

#include <sstream>

template <class Matrix, class Vector>
void solve_PDE_1D(utopia::Input& in) {
    using namespace utopia;

    using Comm = typename Traits<Vector>::Communicator;
    using SizeType = typename Traits<Vector>::SizeType;
    using Scalar = typename Traits<Vector>::Scalar;

    Matrix A;
    Vector x, b;

    // Tipically passed from outside
    Comm& comm = Comm::get_default();

    // Set default values
    SizeType n_local = 10;
    Scalar rhs_value = 1.0;

    // Read from command line
    in.get("n_local", n_local);
    in.get("rhs_value", rhs_value);

    SizeType n_global = n_local * comm.size();
    const Scalar h = 1.0 / n_global;

    auto l = layout(comm, n_local, n_global);

    // Create vectors
    x.zeros(l);
    b.values(l, rhs_value * h);

    // Boundary conditions
    {
        auto r = range(b);
        Write<Vector> w(b);

        if (r.inside(0)) {
            b.set(0, 0.0);
        }

        if (r.inside(n_global - 1)) {
            b.set(n_global - 1, 0.0);
        }
    }

    // Create matrix for 1D FEM Laplacian
    A.sparse(square_matrix_layout(l), 3, 2);

    {
        Write<Matrix> w(A);
        Range r = row_range(A);

        for (SizeType i = r.begin(); i != r.end(); ++i) {
            if (i == 0 || i == n_local - 1) {
                A.set(i, i, 1.0);
                continue;
            }

            if (i > 0) {
                A.set(i, i - 1, -1.0);
            }

            if (i < n_global - 1) {
                A.set(i, i + 1, -1.0);
            }

            A.set(i, i, 2.0);
        }
    }

    A *= h;

    // Solve system
    solve(A, b, x);

    // Compute norm of residual of linear system (with one additional vector allocation)
    Vector r = A * x;
    r = b - r;
    Scalar norm_r = norm2(r);

    std::stringstream ss;
    ss << "norm_r: " << norm_r << std::endl;

    comm.root_print(ss.str());

    // Rename the vector (only supported by few backends)
    rename("x", x);

    // Write to disk
    write("X.m", x);
}

// Run it with `./examples/example_solve_PDE_1D -n_local 10 -rhs_value 1.0`
int main(int argc, char** argv) {
    using namespace utopia;

#ifdef UTOPIA_WITH_PETSC
    using MatrixT = PetscMatrix;
    using VectorT = PetscVector;
#else
#ifdef WITH_TRILINOS
    using MatrixT = TpetraMatrixd;
    using VectorT = TpetraVectord;
#else
    using MatrixT = BlasMatrixd;
    using VectorT = BlasVectord;
#endif
#endif

    Utopia::Init(argc, argv);

    InputParameters params;
    params.init(argc, argv);

    solve_PDE_1D<MatrixT, VectorT>(params);

    return Utopia::Finalize();
}

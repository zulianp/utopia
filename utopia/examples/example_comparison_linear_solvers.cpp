#include <cstdlib>
#include "utopia.hpp"
#include "utopia_BiCGStab.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_GaussSeidel.hpp"
#include "utopia_Jacobi.hpp"
#include "utopia_PointJacobi.hpp"
#include "utopia_Version.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

const char* file = "./linear_solvers.csv";

template <class Solver, class Matrix, class Vector>
void measure_solver(Matrix& A, Vector& b, const string& solver_name) {
    using namespace utopia;
    using Scalar = typename Traits<Vector>::Scalar;

    Solver s;
    s.atol(1e-10);
    s.rtol(1e-10);
    s.stol(1e-10);
    s.max_it(b.size() * 1000);

    Chrono c;
    Vector x(layout(b), 0);

    auto& comm = b.comm();

    c.start();
    // Setting a conditioner is important with iterative linear solver.
    // The more a problem is bad conditioned, that is the more the difference
    // between the largest and small eigenvalues is significant, the more
    // the solution will deviate from the exact solution.
    // A preconditioner allows to make this difference smaller, so that
    // the solution found with the iterative system will be more
    // precise.
    // You can uncomment it.
    // s.set_preconditioner(make_shared<InvDiagPreconditioner<Matrix, Vector>>());

    Vector r = A * x;
    r = b - r;
    x += r;

    s.solve(A, b, x);
    c.stop();

    r = A * x;
    r = b - r;
    Scalar norm_r = norm2(r);

    if (comm.rank() == 0) {
        // stringstream ss;
        cout << solver_name << endl << "norm_r: " << norm_r << " time: " << c.get_seconds() << endl << endl;

        ofstream myfile;
        myfile.open(file, ios_base::app);
        if (!myfile.good()) {
            std::cerr << "Could not open file for writing" << endl;
        } else {
            myfile << solver_name << "," << norm_r << "," << c.get_seconds() << "\n";
            myfile.close();
        }
    }
}

template <class Matrix, class Vector>
void test_linear_solver() {
    using namespace utopia;

    // ---------------------------------------------------------------
    // The code inside here is only to initialise the vectors
    // and the matrix.

    using Comm = typename Traits<Vector>::Communicator;
    using SizeType = typename Traits<Vector>::SizeType;
    using Scalar = typename Traits<Vector>::Scalar;

    Matrix A;
    Vector b, x;

    Comm& comm = Comm::get_default();

    SizeType n_local = 100;
    SizeType n_global = n_local * comm.size();

    auto l = layout(comm, n_local, n_global);

    // Vectors
    b.zeros(l);
    x.zeros(l);

    auto n = b.size();

    {
        auto b_view = local_view_device(b);
        parallel_for(
            local_range_device(b), UTOPIA_LAMBDA(const SizeType& i) {
                const auto bi = b_view.get(i);
                auto h = 1. / Scalar(n);
                auto x = i * h;
                auto f = device::abs(x - 0.5) * h;

                b_view.set(i, f);
            });
    }

    // Matrix
    SizeType n_local_rows = 10;
    SizeType n_local_cols = 10;
    SizeType d_nnz = 3;
    SizeType o_nnz = 3;

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

    // Consider the scaling of the derivative
    // A *= n;

    if (comm.rank() == 0) {
        ofstream ofs;
        ofs.open(file, ofstream::out | ofstream::trunc);
        ofs.close();

        measure_solver<ConjugateGradient<Matrix, Vector>>(A, b, "ConjugateGradient");
        measure_solver<GaussSeidel<Matrix, Vector>>(A, b, "GaussSeidel");
        measure_solver<BiCGStab<Matrix, Vector>>(A, b, "BiCGStab");
        measure_solver<Jacobi<Matrix, Vector>>(A, b, "Jacobi");
        measure_solver<PointJacobi<Matrix, Vector>>(A, b, "PointJacobi");
        measure_solver<GMRES<Matrix, Vector>>(A, b, "GMRES");
        measure_solver<SteihaugToint<Matrix, Vector>>(A, b, "SteihaugToint");

        // measure_solver<Factorization<Matrix, Vector>>(A, b, "Factorization");

        if (system("/usr/bin/python2.7 ./../examples/linear_solvers_plot.py ./linear_solvers.csv 1")) {
            cout << system("/usr/bin/python2.7 ./../examples/linear_solvers_plot.py ./linear_solvers.csv 1");
        } else if (system("/usr/bin/python3 ./../examples/linear_solvers_plot.py ./linear_solvers.csv 1")) {
            cout << system("/usr/bin/python3 ./../examples/linear_solvers_plot.py ./linear_solvers.csv 1");
        } else {
            cout << "Update python please." << endl;
        }
    }
}

// -------------------------------------------------------------

// Run it with `./examples/example_comparison_linear_solvers
int main(const int argc, char* argv[]) {
    using namespace utopia;

// Checking which backend is active.
#ifdef WITH_PETSC
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

    test_linear_solver<MatrixT, VectorT>();

    return Utopia::Finalize();
}

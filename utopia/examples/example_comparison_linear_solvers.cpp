#include "utopia.hpp"
#include "utopia_BiCGStab.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_GaussSeidel.hpp"
#include "utopia_Jacobi.hpp"
#include "utopia_PointJacobi.hpp"
#include "utopia_Version.hpp"
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

struct MySolver {
    string name;
    utopia::Chrono time;
};

// struct less_than_key {
//     inline bool operator() (const MySolver& struct1, const MySolver& struct2) {
//         return (struct1.time.get_seconds() < struct2.time.get_seconds());
//     }
// };

template <class Solver, class Matrix, class Vector>
void measure_solver(Matrix& A, Vector& b, const string solver_name) {
    using namespace utopia;
    using Scalar = typename Traits<Vector>::Scalar;

    Solver s;
    Chrono c;
    Vector x(layout(b), 0);

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
    s.solve(A, b, x);
    c.stop();

    Vector r = A * x;
    r = b - r;
    Scalar norm_r = norm2(r);

    // stringstream ss;
    cout << solver_name << endl << "norm_r: " << norm_r << " time: " << c.get_seconds() << endl << endl;

    std::ofstream myfile;
    myfile.open("./../examples/linear_solvers.csv", ios_base::app);
    if (!myfile.good()) {
        std::cerr << "Could not open file for writing" << std::endl;
    } else {
        myfile << solver_name << ";"<< norm_r<< ";"<< c.get_seconds() << "\n";
        myfile.close();
    }

    


}

template <class Matrix, class Vector>
void test_linear_solver() {
    using namespace utopia;
    vector<MySolver> solvers;

    // ---------------------------------------------------------------
    // The code inside here is only to initialise the vectors
    // and the matrix.

    using Comm = typename Traits<Vector>::Communicator;
    using SizeType = typename Traits<Vector>::SizeType;

    Matrix A;
    Vector b, x;

    Comm& comm = Comm::get_default();

    SizeType n_local = 10;
    SizeType n_global = n_local * comm.size();

    auto l = layout(comm, n_local, n_global);

    // Vectors
    b.zeros(l);
    x.zeros(l);

    {
        auto b_view = local_view_device(b);
        parallel_for(
            local_range_device(b), UTOPIA_LAMBDA(const SizeType& i) {
                const auto bi = b_view.get(i);
                b_view.set(i, i);
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

    // ofstream ofs;
    // ofs.open("./../examples/linear_solver.csv", ofstream::out | ofstream::trunc);
    // ofs.close();

    measure_solver<ConjugateGradient<Matrix, Vector>>(A, b, "ConjugateGradient");
    measure_solver<GaussSeidel<Matrix, Vector>>(A, b, "GaussSeidel");
    measure_solver<BiCGStab<Matrix, Vector>>(A, b, "BiCGStab");
    measure_solver<Jacobi<Matrix, Vector>>(A, b, "Jacobi");
    measure_solver<PointJacobi<Matrix, Vector>>(A, b, "PointJacobi");

    cout << system("/usr/bin/python2.7 ./../examples/hello_world.py 1");

    // // -------------- PRINTING THE RESULTS ---------------------------

    //     sort(solvers.begin(), solvers.end(), less_than_key());
    //     for (auto it = solvers.begin(); it != solvers.end(); ++it){
    //         cout << it->name << ":" << it->time <<endl;
    //     }

    // measure_solver<ConjugateGradient, MatrixT, VectorT>(Matrix &A, Vector &b, "Conjugate Gradient");


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
    // print_statistics_test();

    return Utopia::Finalize();
}

#include "utopia.hpp"
#include "utopia_Version.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;


struct MySolver {
 string name;
 utopia::Chrono time;
};

struct less_than_key {
    inline bool operator() (const MySolver& struct1, const MySolver& struct2) {
        return (struct1.name < struct2.name);
    }
};

template <class Matrix, class Vector>
void fastest_solver() {
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

// End of the initialisation. 
// ---------------------------------------------------------------    
// -------------------- CONJUGATE GRADIENT -----------------------

    Chrono timeCG;
    timeCG.start();
    {
        ConjugateGradient<Matrix, Vector> cg;
        // Setting a conditioner is important with iterative linear solver.
        // The more a problem is bad conditioned, that is the more the difference
        // between the largest and small eigenvalues is significant, the more
        // the solution will deviate from the exact solution. 
        // A preconditioner allows to make this difference smaller, so that
        // the solution found with the iterative system will be more 
        // precise. 
        // Set the preconditioner by uncommenting the following line. 
        //cg.set_preconditioner(make_shared<InvDiagPreconditioner<Matrix, Vector>>());
        cg.solve(A, b, x);
    }
    timeCG.stop();

    solvers.push_back(MySolver());
    solvers[0].name = "Conjugate Gradient";
    solvers[0].time = timeCG;
  
// ---------------------------------------------------------------    
// -------------------- GAUSS SIEDEL -----------------------------

    Chrono timeGS;
    timeGS.start();
    {
        GaussSeidel<Matrix, Vector> gs;
    //  Set the preconditioner by uncommenting the following line. 
    //  (Somehow not working in this example)  
    //  gs.set_preconditioner(make_shared<InvDiagPreconditioner<Matrix, Vector>>());
        gs.solve(A, b, x);
    }
    timeGS.stop();

    solvers.push_back(MySolver());
    solvers[1].name = "Gauss Siedel";
    solvers[1].time = timeGS;

// ---------------------------------------------------------------    
// ------------------------- BiCGStab ------------------------------

    Chrono timeBi;
    timeBi.start();
    {
        BiCGStab<Matrix, Vector> Bi;
        // Set the preconditioner by uncommenting the following line. 
        // Bi.set_preconditioner(make_shared<InvDiagPreconditioner<Matrix, Vector>>());
        Bi.solve(A, b, x);
    }
    timeBi.stop();

    solvers.push_back(MySolver());
    solvers[2].name = "BiCGStab";
    solvers[2].time = timeBi;

// -------------- PRINTING THE RESULTS ---------------------------

    sort(solvers.begin(), solvers.end(), less_than_key());
    for (auto it = solvers.begin(); it != solvers.end(); ++it){
        cout << it->name << ":" << it->time <<endl;
    }

}

// -------------------------------------------------------------

// Run it with `./examples/example_comparison_linear_solvers
int main(const int argc, char *argv[]) {
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

    fastest_solver<MatrixT, VectorT>();

    return Utopia::Finalize();
}


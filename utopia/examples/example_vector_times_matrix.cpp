#include <iostream>
#include <sstream>
#include "utopia.hpp"
// In this simple example we solve the following equation Ab = x.
// A in this case is a simple Identity matrix = 1.
// So in the end we should have that b = x.

template <class Matrix, class Vector>
void vector_by_matrix_mult(utopia::Input& in) {
    using namespace utopia;
    using Comm = typename Traits<Vector>::Communicator;
    using SizeType = typename Traits<Vector>::SizeType;

    // Initialize the Matrix and vectors.
    Matrix A;
    Vector b, x;

    // Setup the communicator.
    Comm& comm = Comm::get_default();

    // Set default values.
    SizeType n_local = 10;
    auto n_value = 2;

    // Give different values to n_value and n_local from command line.
    in.get("n_local", n_local);
    in.get("n_value", n_value);

    // Setup the layout of the vectors and matrices.
    // This allows for everything to communicate properly and set the
    // size of the matrix and vectors to be the same.
    SizeType n_global = n_local * comm.size();
    auto l = layout(comm, n_local, n_global);

    // Set vector b to have a structure as follows: [n_value, n_value,...]*n_local.
    b.values(l, 0);

    {
        auto b_view = local_view_device(b);

        // For loop that iterates through the elements of the vector from index
        // 0 = r.begin() till index r.end() = in this default case 9.
        // Set the vector to have zeros everywhere except at the index 0
        // and index 9.
        // This loop does not take advantage of parallel.
        // for (SizeType i = r.begin(); i != r.end(); ++i) {
        //     if (i > 0 && i < n_local - 1) {
        //         b.set(i, 0);
        //     }
        // }


        // Parallel for loop.
        parallel_for(
            local_range_device(b), UTOPIA_LAMBDA(const SizeType& i) {
                if (i == 0 || i == n_local - 1) {
                    b_view.set(i, n_value);
                }
            });
    }

    A.identity(square_matrix_layout(l), 1);

    // Example of identity matrix.
    // {
    //     // r in a 3X3 matrix has range of 3, the range is how many columns the matrix has. With indexes starting from
    //     0.
    //     // so in this case: 0, 1, 2. Similar to an array
    //     Write<Matrix> w(A);
    //     auto r = row_range(A);

    //     for (SizeType i = r.begin(); i != r.end(); ++i) {
    //         A.set(i, i, 1);
    //     }
    // }

    // Solve the equation.
    solve(A, b, x);

    // Display the vector x.
    disp(x);
}

int main(const int argc, char* argv[]) {
    using namespace utopia;

    // Use certain Matrix or Vector depending
    // on what you have on your machine.
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

    // Initialize utopia...
    Utopia::Init(argc, argv);

    InputParameters params;
    params.init(argc, argv);

    // Display the resulting vector.
    vector_by_matrix_mult<MatrixT, VectorT>(params);

    return Utopia::Finalize();
}
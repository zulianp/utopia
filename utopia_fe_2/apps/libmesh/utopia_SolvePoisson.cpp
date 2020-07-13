#include "utopia_AppRunner.hpp"

#include "utopia.hpp"

#include "utopia_kernels.hpp"

// #include "utopia_petsc_Matrix.hpp"
// #include "utopia_petsc_Vector.hpp"

#include "utopia_libmesh.hpp"

#include <cmath>

namespace utopia {

    static void solve_poisson(Input &in) {
        using Comm = Traits<LMFunctionSpace>::Communicator;
        using Vector = Traits<LMFunctionSpace>::Vector;
        using Matrix = Traits<LMFunctionSpace>::Matrix;
        using Scalar = Traits<LMFunctionSpace>::Scalar;

        Comm comm;
        LMFunctionSpace space(comm);
        space.read(in);

        std::string path = "out.e";
        in.get("output", path);

        Matrix A;
        Vector x, b;

        Assembler<LMFunctionSpace> assembler(space);
        Poisson<Scalar, Scalar> poisson{1.0, 2.0};

        space.create_matrix(A);
        space.create_vector(x);
        space.create_vector(b);

        assembler.assemble(poisson, A, b);
        space.apply_constraints(A, b);

        ConjugateGradient<Matrix, Vector> cg;
        cg.solve(A, b, x);

        space.write(path, x);
    }

    UTOPIA_REGISTER_APP(solve_poisson);

}  // namespace utopia
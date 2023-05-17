// include edsl components
#include "utopia_AppRunner.hpp"

#ifdef UTOPIA_WITH_PETSC
#ifdef UTOPIA_WITH_MATRIX_IO
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

using namespace utopia;

static void convert_matrix_to_raw(Input &args) {
    Path in, out;
    PetscMatrix mat;
    args.require("in", in);
    args.require("out", out);

    mat.read(PetscCommunicator::get_default().get(), in);
    mat.write_raw(PetscCommunicator::get_default().get(), out / "rowptr.raw", out / "colidx.raw", out / "values.raw");
}

UTOPIA_REGISTER_APP(convert_matrix_to_raw);

static void convert_vector_to_raw(Input &args) {
    Path in, out;
    PetscVector vec;
    args.require("in", in);
    args.require("out", out);

    vec.read(PetscCommunicator::get_default().get(), in);
    vec.write_raw(out);
}

UTOPIA_REGISTER_APP(convert_vector_to_raw);

#endif
#endif

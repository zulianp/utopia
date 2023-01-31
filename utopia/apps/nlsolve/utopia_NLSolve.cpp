
#include "utopia_AppRunner.hpp"

#include "utopia_plugin_Function_impl.hpp"

#include "utopia.hpp"

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

void nlsolve(utopia::Input &in) {
    utopia::PluginFunction<Matrix_t> fun;
    fun.read(in);

    utopia::MPICommunicator comm(MPI_COMM_WORLD);
    fun.initialize(comm);
}

UTOPIA_REGISTER_APP(nlsolve);

#include "utopia_petsc_Eval_Monitor.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

namespace utopia {

    void EvalMonitor<PetscMatrix, PETSC>::apply(
        const SizeType &iteration,
        const PetscMatrix &m,
        const std::string &name_of_file,
        const std::string &name_of_mat
        )
    {
        PetscViewer viewer = nullptr;

        PetscViewerASCIIOpen(m.communicator(), name_of_file.c_str(), &viewer);
        PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

        const char *name;
        PetscObjectGetName((PetscObject)m.raw_type(), &name);
        PetscObjectSetName((PetscObject)m.raw_type(), (name_of_mat + std::to_string(iteration)).c_str());
        MatView(m.raw_type(), viewer);

        PetscViewerDestroy(&viewer);
        PetscObjectSetName((PetscObject)m.raw_type(), name);
    }

    void EvalMonitor<PetscVector, PETSC>::apply(
        const SizeType &iteration,
        const PetscVector &v,
        const std::string &name_of_file,
        const std::string &name_of_vec
        )
    {
        PetscViewer viewer = nullptr;

        PetscViewerASCIIOpen(v.communicator(), name_of_file.c_str(), &viewer);
        PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

        const char *name;
        PetscObjectGetName((PetscObject)v.raw_type(), &name);
        PetscObjectSetName((PetscObject)v.raw_type(), (name_of_vec + std::to_string(iteration)).c_str());
        VecView(v.raw_type(), viewer);
        PetscViewerDestroy(&viewer);
        PetscObjectSetName((PetscObject)v.raw_type(), name);
    }

}

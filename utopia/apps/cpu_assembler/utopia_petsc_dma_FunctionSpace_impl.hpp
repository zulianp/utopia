#ifndef UTOPIA_PETSC_DMA_FUNCTIONSPACE_IMPL_HPP
#define UTOPIA_PETSC_DMA_FUNCTIONSPACE_IMPL_HPP

#include "utopia_petsc_dma_FunctionSpace.hpp"
#include "utopia_PetscDM_impl.hpp"

namespace utopia {

    template<class Elem>
    void FunctionSpace<PetscDM<Elem::Dim>, 1, Elem>::elem(const SizeType &idx, Elem &e) const
    {
        mesh_->elem(idx, e);
        typename Mesh::Point translation, cell_size;
        mesh_->cell_point(idx, translation);
        mesh_->cell_size(idx, cell_size);
        e.set(translation, cell_size);
    }

    template<class Elem>
    bool FunctionSpace<PetscDM<Elem::Dim>, 1, Elem>::write(const Path &path, const PetscVector &x) const
    {
        PetscErrorCode ierr = 0;

        PetscViewer       viewer;
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, path.c_str(), &viewer);
        if(ierr != 0) return false;

        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK); assert(ierr == 0);

        // PetscViewerBinaryOpen(PETSC_COMM_WORLD, path.c_str(),FILE_MODE_WRITE, &viewer);
        // PetscViewerPushFormat(viewer, PETSC_VIEWER_BINARY_VTK);

        DMView(raw_type(*mesh_), viewer);
        VecView(raw_type(x), viewer);

        //Extra output
        PetscVector temp = x;
        temp.set(0.0);
        temp.set(x.comm().rank());

        utopia::rename("comm_rank", temp);
        VecView(raw_type(temp), viewer);


        // each_write(temp, [](const SizeType &i) {
        //     return i;
        // });

        // utopia::rename("global_node_id", temp);
        // VecView(raw_type(temp), viewer);

        ierr = PetscViewerDestroy(&viewer); assert(ierr == 0);
        return ierr == 0;
    }

}

#endif //UTOPIA_PETSC_DMA_FUNCTIONSPACE_IMPL_HPP

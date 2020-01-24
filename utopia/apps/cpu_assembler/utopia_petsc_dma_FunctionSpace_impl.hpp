#ifndef UTOPIA_PETSC_DMA_FUNCTIONSPACE_IMPL_HPP
#define UTOPIA_PETSC_DMA_FUNCTIONSPACE_IMPL_HPP

#include "utopia_petsc_dma_FunctionSpace.hpp"
#include "utopia_PetscDM_impl.hpp"
#include "utopia_petsc.hpp"

namespace utopia {

    template<class Elem, int NComponents>
    void FunctionSpace<PetscDM<Elem::Dim>, NComponents, Elem>::elem(const SizeType &idx, Elem &e) const
    {
        mesh_->elem(idx, e.univar_elem());
        typename Mesh::Point translation, cell_size;
        mesh_->cell_point(idx, translation);
        mesh_->cell_size(idx, cell_size);
        e.set(translation, cell_size);
    }

    template<class Elem, int NComponents>
    bool FunctionSpace<PetscDM<Elem::Dim>, NComponents, Elem>::write(const Path &path, const PetscVector &x) const
    {
        PetscErrorCode ierr = 0;
        PetscViewer       viewer;

        const auto ext = path.extension();
        if(ext == "vtk") {
            ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, path.c_str(), &viewer);
            if(ierr != 0) { assert(false); return false; }

            ierr = PetscViewerPushFormat(viewer,  PETSC_VIEWER_ASCII_VTK); assert(ierr == 0);
        } else if(ext == "vts") {

            ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, path.c_str(), FILE_MODE_WRITE, &viewer);
            if(ierr != 0) { assert(false); return false; }

            ierr = PetscViewerPushFormat(viewer,  PETSC_VIEWER_VTK_VTS); assert(ierr == 0);
        } else if(ext == "vtr") {

            ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, path.c_str(), FILE_MODE_WRITE, &viewer);
            if(ierr != 0) { assert(false); return false; }

            ierr = PetscViewerPushFormat(viewer,  PETSC_VIEWER_VTK_VTR); assert(ierr == 0);
        } else if(ext == "vtu") {

            ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, path.c_str(), FILE_MODE_WRITE, &viewer);
            if(ierr != 0) { assert(false); return false; }

            ierr = PetscViewerPushFormat(viewer,  PETSC_VIEWER_VTK_VTU); assert(ierr == 0);
        } else {
            std::cerr << "unknown format " << ext << std::endl;
            return false;
        }

        // if(NComponents == 1 && mesh_->n_components() > 1) {
        //     const SizeType n = x.local_size()/mesh_->n_components();

        //     PetscVector s = local_zeros(n);
        //     ierr = VecStrideGather(raw_type(x), subspace_id_, raw_type(s), INSERT_VALUES); assert(ierr == 0);

        //     ierr = DMView(raw_type(*mesh_), viewer);    assert(ierr == 0);
        //     ierr = VecView(raw_type(s), viewer);        assert(ierr == 0);
        // } else {
            ierr = DMView(raw_type(*mesh_), viewer);    assert(ierr == 0);
            ierr = VecView(raw_type(x), viewer);        assert(ierr == 0);
        // }

        ierr = PetscViewerPopFormat(viewer); assert(ierr == 0);
        ierr = PetscViewerDestroy(&viewer);  assert(ierr == 0);
        return ierr == 0;
    }

}

#endif //UTOPIA_PETSC_DMA_FUNCTIONSPACE_IMPL_HPP

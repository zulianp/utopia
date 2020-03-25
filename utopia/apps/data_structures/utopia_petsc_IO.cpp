
#include "utopia_petsc_IO.hpp"
#include "utopia_Path.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_DM.hpp"

#include <iostream>
#include <petscviewerhdf5.h>

namespace utopia {

    class PetscIO::Wrapper {
    public:
        Wrapper()
        : viewer(nullptr)
        {}

        ~Wrapper()
        {
            destroy();
        }

        void destroy()
        {
            if(viewer) {
                auto ierr = PetscViewerDestroy(&viewer);  assert(ierr == 0);
                viewer = nullptr;
            }
        }

        PetscViewer viewer;
    };

    bool PetscIO::open(const PetscCommunicator &comm, const Path &path)
    {
        PetscErrorCode ierr = 0;

        auto mpi_comm = comm.get();
        auto &viewer = wrapper_->viewer;

        const auto ext = path.extension();
        if(ext == "vtk") {
            ierr = PetscViewerASCIIOpen(mpi_comm, path.c_str(), &viewer);
            if(ierr != 0) { assert(false); return false; }

            ierr = PetscViewerPushFormat(viewer,  PETSC_VIEWER_ASCII_VTK); assert(ierr == 0);
        } else if(ext == "vts") {

            ierr = PetscViewerVTKOpen(mpi_comm, path.c_str(), FILE_MODE_WRITE, &viewer);
            if(ierr != 0) { assert(false); return false; }

            ierr = PetscViewerPushFormat(viewer,  PETSC_VIEWER_VTK_VTS); assert(ierr == 0);
        } else if(ext == "vtr") {

            ierr = PetscViewerVTKOpen(mpi_comm, path.c_str(), FILE_MODE_WRITE, &viewer);
            if(ierr != 0) { assert(false); return false; }

            ierr = PetscViewerPushFormat(viewer,  PETSC_VIEWER_VTK_VTR); assert(ierr == 0);
        } else if(ext == "vtu") {

            ierr = PetscViewerVTKOpen(mpi_comm, path.c_str(), FILE_MODE_WRITE, &viewer);
            if(ierr != 0) { assert(false); return false; }

            ierr = PetscViewerPushFormat(viewer,  PETSC_VIEWER_VTK_VTU); assert(ierr == 0);
        }
#if defined(PETSC_HAVE_HDF5)
        else if(ext == "h5") {
            PetscViewerHDF5Open(mpi_comm, path.c_str(), FILE_MODE_WRITE, &viewer);
            return true;
        }
#endif
        else if(ext == "png") {
            //It would be nice to support direct rendering for 2D
            std::cerr << "unsupported format " << ext << std::endl;
            return false;
        } else {
            std::cerr << "unknown format " << ext << std::endl;
            return false;
        }


        return true;
    }

    bool PetscIO::write(const PetscDMBase &dm)
    {
        return DMView(dm.raw_type(), wrapper_->viewer) == 0;
    }

    bool PetscIO::write(const PetscVector &x)
    {
        return VecView(x.raw_type(), wrapper_->viewer) == 0;
    }

    void PetscIO::close()
    {
        wrapper_->destroy();
    }

    PetscIO::PetscIO() : wrapper_(utopia::make_unique<Wrapper>())
    {}

    PetscIO::~PetscIO() {}

}


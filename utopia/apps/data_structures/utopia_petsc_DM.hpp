#ifndef UTOPIA_PETSC_DM_HPP
#define UTOPIA_PETSC_DM_HPP

#include "utopia_Input.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_petsc_Communicator.hpp"
#include "utopia_Path.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Matrix.hpp"

#include <petscdm.h>
#include <petscviewerhdf5.h>

namespace utopia {
    class PetscDMBase : public Configurable {
    public:
        virtual ~PetscDMBase() {}

        template<class SizeType>
        static void dof_ownership_range(const DM &dm, SizeType &begin, SizeType &end)
        {
            Vec v;
            DMGetGlobalVector(dm, &v);
            VecGetOwnershipRange(v, &begin, &end);
            DMRestoreGlobalVector(dm, &v);
        }

        static PetscInt get_dimension(DM dm)
        {
            PetscInt ret;
            DMGetDimension(dm, &ret);
            return ret;
        }

        static void refine(const DM &in, const MPI_Comm &comm, DM &out)
        {
            DMRefine(in, comm, &out);
        }

        DM &raw_type()
        {
            return wrapper_->dm;
        }

        const DM &raw_type() const
        {
            return wrapper_->dm;
        }

        virtual void wrap(DM &dm, const bool delegate_ownership)
        {
            destroy_dm();
            wrapper_->dm = dm;
            wrapper_->owned = delegate_ownership;
        }

        void destroy_dm()
        {
            wrapper_->destroy();
        }

        PetscDMBase(const PetscCommunicator &comm = PetscCommunicator())
        : comm_(comm), wrapper_(utopia::make_unique<Wrapper>())
        {}

        inline PetscCommunicator &comm()
        {
            return comm_;
        }

        inline const PetscCommunicator &comm() const
        {
            return comm_;
        }

        void create_matrix(PetscMatrix &mat) const
        {
            mat.destroy();
            DMCreateMatrix(raw_type(), &mat.raw_type());
        }

        void create_vector(PetscVector &vec) const
        {
            vec.destroy();
            DMCreateGlobalVector(raw_type(), &vec.raw_type());
        }

        void create_local_vector(PetscVector &vec) const
        {
            vec.destroy();
            auto err = DMCreateLocalVector(raw_type(), &vec.raw_type()); assert(err == 0);
        }

        void create_interpolation(const PetscDMBase &target, PetscMatrix &I) const
        {
            I.destroy();
            auto ierr = DMCreateInterpolation(raw_type(), target.raw_type(), &I.raw_type(), nullptr); assert(ierr == 0);
        }

        void local_to_global(const PetscVector &local,  PetscVector &global) const
        {
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 0) //DM-INCOMPLETE
            DMLocalToGlobal(raw_type(), local.raw_type(), ADD_VALUES, global.raw_type());
#else
            DMLocalToGlobalBegin(raw_type(), local.raw_type(), ADD_VALUES, global.raw_type());
            DMLocalToGlobalEnd(raw_type(), local.raw_type(), ADD_VALUES, global.raw_type());
#endif
        }

        void global_to_local(const PetscVector &global, PetscVector &local) const
        {
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 0) //DM-INCOMPLETE
            DMGlobalToLocal(raw_type(), global.raw_type(), INSERT_VALUES, local.raw_type());
#else
            DMGlobalToLocalBegin(raw_type(), global.raw_type(), INSERT_VALUES, local.raw_type());
            DMGlobalToLocalEnd(raw_type(), global.raw_type(), INSERT_VALUES, local.raw_type());
#endif
        }

        bool write(const Path &path, const PetscVector &x) const
        {
            PetscErrorCode ierr = 0;
            PetscViewer    viewer;

            auto mpi_comm = comm().get();

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

                ierr = DMView(raw_type(*mesh_), viewer);    assert(ierr == 0);
                ierr = VecView(raw_type(x), viewer);        assert(ierr == 0);

                ierr = PetscViewerDestroy(&viewer);  assert(ierr == 0);
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


            ierr = DMView(this->raw_type(), viewer);    assert(ierr == 0);
            ierr = VecView(x.raw_type(), viewer);        assert(ierr == 0);

            ierr = PetscViewerPopFormat(viewer); assert(ierr == 0);
            ierr = PetscViewerDestroy(&viewer);  assert(ierr == 0);
            return ierr == 0;
        }

    private:
        PetscCommunicator comm_;

        class Wrapper {
        public:
            Wrapper()
            : dm(nullptr), owned(true)
            {}

            ~Wrapper() {
                destroy();
            }

            void destroy()
            {
                if(owned && dm) {
                    DMDestroy(&dm);
                    dm = nullptr;
                }
            }

            DM dm;
            bool owned;
        };

        std::unique_ptr<Wrapper> wrapper_;
    };
}

#endif //UTOPIA_PETSC_DM_HPP

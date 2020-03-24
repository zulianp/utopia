#ifndef UTOPIA_PETSC_DMDA_HPP
#define UTOPIA_PETSC_DMDA_HPP

#include "utopia_StructuredGrid.hpp"
#include "utopia_petsc_DM.hpp"

#include <petscdmda.h>
#include <cassert>

namespace utopia {

    template<class Point, class IntArray>
    class PetscDMDA : public StructuredGrid<Point, IntArray>, public PetscDMBase {
    public:
        using Super = utopia::StructuredGrid<Point, IntArray>;

        using SizeType = typename Super::SizeType;
        using Scalar   = typename Super::Scalar;

        void init(DM dm)
        {
            PetscInt dof_range_begin, dof_range_end;
            PetscDMBase::dof_ownership_range(dm, dof_range_begin, dof_range_end);

            const SizeType dim = PetscDMBase::get_dimension(dm);
            assert(dim == this->dim());

            this->set_dof_range_begin(dof_range_begin);
            this->set_dof_range_end(dof_range_end);

            get_dims(dm, this->dims());
            get_corners(dm, this->corners_begin(), this->corners_extent());
            get_ghost_corners(dm, this->ghost_corners_begin(), this->ghost_corners_extent());
            this->set_n_components( get_dof(dm) );

            DMDAElementType elem_type;
            DMDAGetElementType(dm, &elem_type);

            if(elem_type == DMDA_ELEMENT_P1) {
                if(this->dim() == 2) {
                    this->set_elements_x_cell(2);
                } else if(this->dim() == 3) {
                    this->set_elements_x_cell(6);
                } else {
                    assert(false);
                }
            }
        }

        PetscDMDA()
        {
            this->box_min().set(0.0);
            this->box_max().set(1.0);
            this->set_elements_x_cell(1);
        }

        ////////////////////////////////////////////////////////

        static void get_dims(DM dm, IntArray &dims)
        {
            PetscErrorCode ierr;

            utopia::ArrayView<PetscInt, 3> dims_buff;
            ierr = DMDAGetInfo(dm, nullptr, &dims_buff[0], &dims_buff[1], &dims_buff[2], nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); assert(ierr == 0);

            const SizeType n = dims.size();
            for(SizeType d = 0; d < n; ++d) {
                dims[d] = dims_buff[d];
            }
        }

        static void get_corners(DM dm, IntArray &start, IntArray &extent)
        {
            PetscErrorCode ierr;

            utopia::ArrayView<PetscInt, 3> start_buff;
            utopia::ArrayView<PetscInt, 3> extent_buff;

            ierr = DMDAGetCorners(dm,
                &start_buff[0],  &start_buff[1],  &start_buff[2],
                &extent_buff[0], &extent_buff[1], &extent_buff[2]
            ); assert(ierr == 0);

            const SizeType n = start.size();
            for(SizeType d = 0; d < n; ++d) {
                start[d]  = start_buff[d];
                extent[d] = extent_buff[d];
            }
        }

        static void get_ghost_corners(DM dm, IntArray &start, IntArray &extent)
        {
            PetscErrorCode ierr;

            utopia::ArrayView<PetscInt, 3> start_buff;
            utopia::ArrayView<PetscInt, 3> extent_buff;

            ierr = DMDAGetGhostCorners(dm,
                &start_buff[0],  &start_buff[1],  &start_buff[2],
                &extent_buff[0], &extent_buff[1], &extent_buff[2]
            ); assert(ierr == 0);

            const SizeType n = start.size();

            for(SizeType d = 0; d < n; ++d) {
                start[d]  = start_buff[d];
                extent[d] = extent_buff[d];
            }
        }

        static SizeType get_dof(DM dm)
        {
            PetscErrorCode ierr;

            SizeType ret;
            ierr = DMDAGetDof(dm, &ret); assert(ierr == 0);
            return ret;
        }

    };

}

#endif //UTOPIA_PETSC_DMDA_HPP

#ifndef UTOPIA_PETSC_DMDA_HPP
#define UTOPIA_PETSC_DMDA_HPP

#include "utopia_StructuredGrid.hpp"
#include "utopia_petsc_DM.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_Algorithms.hpp"

#include <cassert>

#include <petscdmda.h>


namespace utopia {

    template<class Point, class IntArray>
    class PetscDMDA;

    template<class Point, class IntArray>
    class PetscDMDA : public StructuredGrid<Point, IntArray>, public PetscDMBase {
    public:
        using Super = utopia::StructuredGrid<Point, IntArray>;

        using SizeType  = typename Super::SizeType;
        using Scalar    = typename Super::Scalar;
        using NodeIndex = utopia::ArrayView<SizeType>;

        class Elements {
        public:
            Elements(DM dm)
            : dm_(dm)
            {
                DMDAGetElements(dm_, &n_local_elem_, &n_nodes_x_elem_, &local_elem_);
            }

            ~Elements()
            {
                DMDARestoreElements(dm_, &n_local_elem_, &n_nodes_x_elem_, &local_elem_);
            }

            inline SizeType local_size() const
            {
                return n_local_elem_;
            }

            inline NodeIndex nodes_local(const SizeType &local_elem_idx) const
            {
                return NodeIndex(
                    &local_elem_[local_elem_idx * n_nodes_x_elem_],
                    n_nodes_x_elem_
                );
            }

        private:
            DM dm_;
            SizeType n_local_elem_, n_nodes_x_elem_;
            const SizeType *local_elem_;
        };

        std::unique_ptr<Elements> make_elements() const
        {
            return utopia::make_unique<Elements>(this->raw_type());
        }

        void wrap(DM &dm, const bool delegate_ownership) override
        {
            PetscDMBase::wrap(dm, delegate_ownership);
            init_from_dm(dm);
        }

        PetscDMDA(const PetscCommunicator &comm, const DMDAElementType &type_override = DMDA_ELEMENT_Q1)
        : PetscDMBase(comm), type_override_(type_override) {}

        PetscDMDA(DM &dm, const bool delegate_ownership)
        {
            wrap(dm, delegate_ownership);
        }

        void read(Input &in) override
        {
            init_default();

            const SizeType n = this->dims().size();
            assert(n > 0); //IMPLEMENT ME for dynamically sized arrays

            const char coord_names[3] = {'x', 'y', 'z'};

            std::string n_str   = "n?";
            std::string min_str = "?_min";
            std::string max_str = "?_max";

            for(SizeType d = 0; d < n; ++d) {
                const char coord = coord_names[d];

                n_str[1]   = coord;
                min_str[0] = coord;
                max_str[0] = coord;

                in.get(n_str,   this->dims()[d]);
                in.get(min_str, this->box_min()[d]);
                in.get(max_str, this->box_max()[d]);
            }

            this->destroy_dm();
            create_uniform(
                comm().get(),
                this->dims(),
                this->box_min(),
                this->box_max(),
                type_override_,
                this->n_components(),
                this->raw_type()
            );

            ///re-initialize mirror just to be safe (it is cheap anyway)
            update_mirror();
        }

        void update_mirror()
        {
            init_from_dm(raw_type());
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        static void create(
            const MPI_Comm &comm,
            const IntArray &dims,
            const SizeType n_components,
            DMDAStencilType stencil_type,
            DM &dm)
        {
            const SizeType dim = dims.size();

            switch(dim)
            {
                case 1:
                {
                    DMDACreate1d(comm,
                        DM_BOUNDARY_NONE,
                        dims[0],
                        n_components, //dofs per node
                        1, //stencil width
                        nullptr,
                        &dm
                    );

                    break;
                }

                case 2:
                {
                    DMDACreate2d(comm,
                        DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
                        stencil_type,
                        dims[0], dims[1],
                        PETSC_DECIDE, PETSC_DECIDE,
                        n_components, //dofs per node
                        1, //stencil width
                        nullptr,
                        nullptr,
                        &dm
                    );

                    break;
                }

                case 3:
                {
                    DMDACreate3d(comm,
                        DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
                        stencil_type,
                        dims[0], dims[1], dims[2],
                        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                        n_components, //dofs per node
                        1, //stencil width
                        nullptr,
                        nullptr,
                        nullptr,
                        &dm
                    );

                    break;
                }

                default:
                {
                    assert(false);
                }
            }
        }

        static void create_uniform(
            const MPI_Comm &comm,
            const IntArray &dims,
            const Point &box_min,
            const Point &box_max,
            DMDAElementType elem_type,
            const SizeType n_components,
            DM &dm,
            DMDAStencilType stencil_type = DMDA_STENCIL_BOX
        )
        {
            create(comm, dims, n_components, stencil_type, dm);

            Scalar min_x = 0, min_y = 0, min_z = 0;
            Scalar max_x = 1, max_y = 1, max_z = 1;

            min_x = box_min[0];
            max_x = box_max[0];

            const SizeType d = dims.size();

            if(d > 1) {
                min_y = box_min[1];
                max_y = box_max[1];
            }

            if(d > 2) {
                min_z = box_min[2];
                max_z = box_max[2];
            }

            DMDASetElementType(dm, elem_type);
            DMSetUp(dm);
            DMDASetUniformCoordinates(dm, min_x, max_x, min_y, max_y, min_z, max_z);

            if(elem_type == DMDA_ELEMENT_Q1) {
                DMDASetInterpolationType(dm, DMDA_Q1);
            }
        }

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

        std::unique_ptr<PetscDMDA> uniform_refine() const
        {
            auto fine = utopia::make_unique<PetscDMDA>();
            PetscDMBase::refine(raw_type(), comm().get(), fine->raw_type());

            //This does not transfer automatically for some reason
            DMDAElementType elem_type;
            DMDAGetElementType(raw_type(), &elem_type);
            DMDASetElementType(fine->raw_type(), elem_type);
            fine->update_mirror();

            return std::move(fine);
        }

        std::unique_ptr<PetscDMDA> clone(const SizeType &n_components) const
        {
            auto cloned = utopia::make_unique<PetscDMDA>();
            cloned->copy(*this);
            cloned->set_n_components(n_components);
            cloned->type_override_ = type_override_;
            cloned->init_from_mirror();
            return std::move(cloned);
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

    private:
        DMDAElementType type_override_;

        void init_default()
        {
            device::fill(10, this->dims());
            device::fill(0,  this->box_min());
            device::fill(1,  this->box_max());
            this->set_n_components(1);
        }

        void init_from_mirror()
        {
            this->destroy_dm();
            create_uniform(
                comm().get(),
                this->dims(),
                this->box_min(),
                this->box_max(),
                type_override_,
                this->n_components(),
                this->raw_type()
            );
        }

        void init_from_dm(DM dm)
        {
            MPI_Comm mpi_comm = PetscObjectComm((PetscObject) dm);
            comm().set(mpi_comm);

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

    };

}

#endif //UTOPIA_PETSC_DMDA_HPP

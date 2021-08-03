#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMBase.hpp"

#include "utopia_mesh_StructuredGrid.hpp"

#include <petscdmda.h>

namespace utopia {

    // template <typename T>
    // class Traits<std::vector<T>> {
    // public:
    //     using Scalar = T;
    // };

    namespace petsc {

        class StructuredGrid::Impl : public DMBase,
                                     public mesh::StructuredGrid<std::vector<PetscScalar>, std::vector<PetscInt>> {
        public:
            static void get_corners(DM dm, IntArray &start, IntArray &extent) {
                PetscErrorCode ierr;

                utopia::ArrayView<PetscInt, 3> start_buff{};
                utopia::ArrayView<PetscInt, 3> extent_buff{};

                ierr = DMDAGetCorners(dm,
                                      &start_buff[0],
                                      &start_buff[1],
                                      &start_buff[2],
                                      &extent_buff[0],
                                      &extent_buff[1],
                                      &extent_buff[2]);
                assert(ierr == 0);

                const SizeType n = start.size();
                for (SizeType d = 0; d < n; ++d) {
                    start[d] = start_buff[d];
                    extent[d] = extent_buff[d];
                }

                UTOPIA_UNUSED(ierr);
            }

            static void get_ghost_corners(DM dm, IntArray &start, IntArray &extent) {
                PetscErrorCode ierr;

                utopia::ArrayView<PetscInt, 3> start_buff{};
                utopia::ArrayView<PetscInt, 3> extent_buff{};

                ierr = DMDAGetGhostCorners(dm,
                                           &start_buff[0],
                                           &start_buff[1],
                                           &start_buff[2],
                                           &extent_buff[0],
                                           &extent_buff[1],
                                           &extent_buff[2]);
                assert(ierr == 0);

                const SizeType n = start.size();

                for (SizeType d = 0; d < n; ++d) {
                    start[d] = start_buff[d];
                    extent[d] = extent_buff[d];
                }

                UTOPIA_UNUSED(ierr);
            }

            static SizeType get_dof(DM dm) {
                PetscErrorCode ierr;

                SizeType ret;
                ierr = DMDAGetDof(dm, &ret);
                assert(ierr == 0);
                UTOPIA_UNUSED(ierr);
                return ret;
            }

            class Elements {
            public:
                Elements(DM dm) : dm_(dm) {
                    DMDAGetElements(dm_, &n_local_elem_, &n_nodes_x_elem_, &local_elem_);

                    const SizeType n_idx = n_local_elem_ * n_nodes_x_elem_;

                    global_elem_.resize(n_idx);

                    ISLocalToGlobalMapping map;
                    DMGetLocalToGlobalMapping(dm, &map);
                    ISLocalToGlobalMappingApplyBlock(map, n_idx, local_elem_, &global_elem_[0]);
                }

                ~Elements() { DMDARestoreElements(dm_, &n_local_elem_, &n_nodes_x_elem_, &local_elem_); }

                inline SizeType local_size() const { return n_local_elem_; }
                inline Range range() const { return Range(0, n_local_elem_); }

                inline NodeIndex nodes(const SizeType &local_elem_idx) const {
                    return NodeIndex(&global_elem_[local_elem_idx * n_nodes_x_elem_], n_nodes_x_elem_);
                }

                inline NodeIndex nodes_local(const SizeType &local_elem_idx) const {
                    return NodeIndex(&local_elem_[local_elem_idx * n_nodes_x_elem_], n_nodes_x_elem_);
                }

                inline SizeType first_node_idx(const SizeType &local_elem_idx) const {
                    return local_elem_[local_elem_idx * n_nodes_x_elem_];
                }

            private:
                DM dm_;
                SizeType n_local_elem_, n_nodes_x_elem_;
                const SizeType *local_elem_;
                std::vector<SizeType> global_elem_;
            };

            void init_elements() {
                if (!elements_) {
                    elements_ = make_elements();
                }
            }

            std::shared_ptr<Elements> elements_ptr() {
                init_elements();
                return elements_;
            }

            std::unique_ptr<Elements> make_elements() const { return utopia::make_unique<Elements>(this->raw_type()); }

            static void create(const MPI_Comm &comm,
                               const IntArray &dims,
                               const SizeType n_components,
                               DMDAStencilType stencil_type,
                               DM &dm) {
                const SizeType dim = dims.size();

                switch (dim) {
                    case 1: {
                        DMDACreate1d(comm,
                                     DM_BOUNDARY_NONE,
                                     dims[0],
                                     n_components,  // dofs per node
                                     1,             // stencil width
                                     nullptr,
                                     &dm);

                        break;
                    }

                    case 2: {
                        DMDACreate2d(comm,
                                     DM_BOUNDARY_NONE,
                                     DM_BOUNDARY_NONE,
                                     stencil_type,
                                     dims[0],
                                     dims[1],
                                     PETSC_DECIDE,
                                     PETSC_DECIDE,
                                     n_components,  // dofs per node
                                     1,             // stencil width
                                     nullptr,
                                     nullptr,
                                     &dm);

                        break;
                    }

                    case 3: {
                        DMDACreate3d(comm,
                                     DM_BOUNDARY_NONE,
                                     DM_BOUNDARY_NONE,
                                     DM_BOUNDARY_NONE,
                                     stencil_type,
                                     dims[0],
                                     dims[1],
                                     dims[2],
                                     PETSC_DECIDE,
                                     PETSC_DECIDE,
                                     PETSC_DECIDE,
                                     n_components,  // dofs per node
                                     1,             // stencil width
                                     nullptr,
                                     nullptr,
                                     nullptr,
                                     &dm);

                        break;
                    }

                    default: {
                        assert(false);
                    }
                }
            }

            static void create_uniform(const MPI_Comm &comm,
                                       const IntArray &dims,
                                       const Point &box_min,
                                       const Point &box_max,
                                       DMDAElementType elem_type,
                                       const SizeType n_components,
                                       DM &dm,
                                       DMDAStencilType stencil_type = DMDA_STENCIL_BOX) {
                create(comm, dims, n_components, stencil_type, dm);

                Scalar min_x = 0, min_y = 0, min_z = 0;
                Scalar max_x = 1, max_y = 1, max_z = 1;

                min_x = box_min[0];
                max_x = box_max[0];

                const SizeType d = dims.size();

                if (d > 1) {
                    min_y = box_min[1];
                    max_y = box_max[1];
                }

                if (d > 2) {
                    min_z = box_min[2];
                    max_z = box_max[2];
                }

                DMDASetElementType(dm, elem_type);
                DMSetUp(dm);
                DMDASetUniformCoordinates(dm, min_x, max_x, min_y, max_y, min_z, max_z);

                if (elem_type == DMDA_ELEMENT_Q1) {
                    DMDASetInterpolationType(dm, DMDA_Q1);
                }
            }

            static void get_dims(DM dm, IntArray &dims) {
                PetscErrorCode ierr;

                utopia::ArrayView<PetscInt, 3> dims_buff{};
                ierr = DMDAGetInfo(dm,
                                   nullptr,
                                   &dims_buff[0],
                                   &dims_buff[1],
                                   &dims_buff[2],
                                   nullptr,
                                   nullptr,
                                   nullptr,
                                   nullptr,
                                   nullptr,
                                   nullptr,
                                   nullptr,
                                   nullptr,
                                   nullptr);
                assert(ierr == 0);

                const SizeType n = dims.size();
                for (SizeType d = 0; d < n; ++d) {
                    dims[d] = dims_buff[d];
                }

                UTOPIA_UNUSED(ierr);
            }

            void init_default() {
                elements_ = nullptr;
                device::fill(10, this->dims());
                device::fill(0, this->box_min());
                device::fill(1, this->box_max());

                init_elem_x_cell(type_override_);
            }

            void init_from_mirror() {
                this->destroy_dm();

                create_uniform(comm().get(),
                               this->dims(),
                               this->box_min(),
                               this->box_max(),
                               type_override_,
                               this->n_components(),
                               this->raw_type());

                init_elem_x_cell(type_override_);
            }

            void init_from_dm(DM dm) {
                elements_ = nullptr;

                MPI_Comm mpi_comm = PetscObjectComm((PetscObject)dm);
                comm().set(mpi_comm);

                PetscInt dof_range_begin, dof_range_end;
                DMBase::dof_ownership_range(dm, dof_range_begin, dof_range_end);

#ifndef NDEBUG
                {
                    const SizeType dim = DMBase::get_dimension(dm);
                    assert(dim == this->dim());
                }
#endif  // NDEBUG

                this->set_dof_range_begin(dof_range_begin);
                this->set_dof_range_end(dof_range_end);

                get_dims(dm, this->dims());
                get_corners(dm, this->corners_begin(), this->corners_extent());
                get_ghost_corners(dm, this->ghost_corners_begin(), this->ghost_corners_extent());
                this->set_n_components(get_dof(dm));

                DMDAElementType elem_type;
                DMDAGetElementType(dm, &elem_type);
                init_elem_x_cell(elem_type);
            }

            void init_elem_x_cell(DMDAElementType elem_type) {
                int exc = 1;
                if (elem_type == DMDA_ELEMENT_P1) {
                    if (this->dim() == 2) {
                        exc = 2;
                    } else if (this->dim() == 3) {
                        exc = 6;
                    } else {
                        assert(false);
                    }
                }

                this->set_elements_x_cell(exc);
            }

            // IntArray &dims() { return dims_; }
            // Point &box_min() {}

            // Point &box_max() {}

            // SizeType n_components() {}

            // inline void cell_point(const SizeType &idx, Point &p) const {
            //     assert(elements_ && "init_elements() has to be called before any access to element data");
            //     const SizeType p_idx = elements_->first_node_idx(idx);
            //     this->point(p_idx, p);
            // }

            void wrap(DM &dm, const bool delegate_ownership) {
                DMBase::wrap(dm, delegate_ownership);
                init_from_dm(dm);
            }

            Impl(const PetscCommunicator &comm, const DMDAElementType &type_override = DMDA_ELEMENT_Q1)
                : DMBase(comm), type_override_(type_override) {
                this->set_n_components(1);
            }

            Impl(DM &dm, const bool delegate_ownership) {
                this->set_n_components(1);
                wrap(dm, delegate_ownership);
            }

            void read(Input &in) override { assert(false); }

            DMDAElementType type_override_;

            // initialized in a lazy way
            std::shared_ptr<Elements> elements_;
            IntArray dims_;
        };

        StructuredGrid::Communicator &StructuredGrid::comm() { return impl_->comm(); }
        const StructuredGrid::Communicator &StructuredGrid::comm() const { return impl_->comm(); }

        std::unique_ptr<StructuredGrid> StructuredGrid::uniform_refine() const {
            auto fine = utopia::make_unique<StructuredGrid>();
            fine->impl_ = utopia::make_unique<StructuredGrid::Impl>(comm(), impl_->type_override_);
            DMBase::refine(impl_->raw_type(), comm().get(), fine->impl_->raw_type());

            device::copy(this->impl_->box_min(), fine->impl_->box_min());
            device::copy(this->impl_->box_max(), fine->impl_->box_max());

            // This does not transfer automatically for some reason
            DMDAElementType elem_type;
            DMDAGetElementType(impl_->raw_type(), &elem_type);
            DMDASetElementType(fine->impl_->raw_type(), elem_type);
            fine->update_mirror();

            return fine;
        }

        std::unique_ptr<StructuredGrid> StructuredGrid::clone(const SizeType &n_components) const {
            auto cloned = utopia::make_unique<StructuredGrid>();
            cloned->impl_ = utopia::make_unique<StructuredGrid::Impl>(comm(), impl_->type_override_);
            cloned->impl_->copy(*this->impl_);
            cloned->impl_->set_n_components(n_components);
            // cloned->type_override_ = type_override_;
            cloned->impl_->init_from_mirror();
            return cloned;
        }

        std::unique_ptr<StructuredGrid> StructuredGrid::clone() const { return clone(this->impl_->n_components()); }

        void StructuredGrid::read(Input &in) {
            impl_->init_default();

            const SizeType n = this->impl_->dims().size();
            assert(n > 0);  // IMPLEMENT ME for dynamically sized arrays

            const char coord_names[3] = {'x', 'y', 'z'};

            std::string n_str = "n?";
            std::string min_str = "?_min";
            std::string max_str = "?_max";

            for (SizeType d = 0; d < n; ++d) {
                const char coord = coord_names[d];

                n_str[1] = coord;
                min_str[0] = coord;
                max_str[0] = coord;

                in.get(n_str, this->impl_->dims()[d]);
                in.get(min_str, this->impl_->box_min()[d]);
                in.get(max_str, this->impl_->box_max()[d]);
            }

            this->impl_->destroy_dm();
            impl_->create_uniform(comm().get(),
                                  this->impl_->dims(),
                                  this->impl_->box_min(),
                                  this->impl_->box_max(),
                                  this->impl_->type_override_,
                                  this->impl_->n_components(),
                                  this->impl_->raw_type());

            /// re-initialize mirror just to be safe (it is cheap anyway)
            update_mirror();
        }

        void StructuredGrid::set_field_name(const SizeType &nf, const std::string &name) {
            DMDASetFieldName(impl_->raw_type(), nf, name.c_str());
        }

        void StructuredGrid::set_field_names(const std::vector<std::string> &names) {
            const std::size_t n_names = names.size();
            std::vector<const char *> names_copy(n_names + 1);

            for (std::size_t i = 0; i < n_names; ++i) {
                names_copy[i] = names[i].c_str();
            }

            names_copy[n_names] = nullptr;

            DMDASetFieldNames(impl_->raw_type(), &names_copy[0]);
        }

        void StructuredGrid::nodes(const SizeType &idx, NodeIndex &nodes) const {
            assert(elements_ && "init_elements() has to be called before any access to element data");
            nodes = impl_->elements_->nodes(idx);
        }

        void StructuredGrid::nodes_local(const SizeType &idx, NodeIndex &nodes) const {
            assert(elements_ && "init_elements() has to be called before any access to element data");
            nodes = impl_->elements_->nodes_local(idx);
        }

        bool StructuredGrid::on_boundary(const SizeType &elem_idx) const {
            NodeIndex idx;
            nodes(elem_idx, idx);

            for (auto i : idx) {
                if (this->impl_->is_node_on_boundary(i)) {
                    return true;
                }
            }

            return false;
        }

        void StructuredGrid::update_mirror() { impl_->init_from_dm(impl_->raw_type()); }

        void StructuredGrid::describe(std::ostream &os) const {
            os << "n_elements      : " << this->impl_->n_elements() << '\n';
            os << "n_nodes         : " << this->impl_->n_nodes() << '\n';
            os << "dim             : " << this->impl_->dim() << '\n';
            os << "elements_x_cell : " << this->impl_->elements_x_cell() << '\n';

            // disp("box_min");
            // disp(this->box_min());

            // disp("box_max");
            // disp(this->box_max());
        }

    }  // namespace petsc
}  // namespace utopia
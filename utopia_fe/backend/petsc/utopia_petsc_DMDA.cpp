#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMBase.hpp"

#include "utopia_mesh_StructuredGrid.hpp"

#include <petscdmda.h>

namespace utopia {
    namespace petsc {

        class DMDABase : public DMBase {
        public:
            using Traits = utopia::Traits<PetscVector>;
            using Communicator = Traits::Communicator;
            using SizeType = Traits::SizeType;
            using Scalar = Traits::Scalar;
            using NodeIndex = utopia::ArrayView<const SizeType>;
            using IntArray = utopia::ArrayView<SizeType, 3>;
            using Point = utopia::StaticVector<Scalar, 3>;

            DMDABase() = default;
            DMDABase(const PetscCommunicator &comm) : DMBase(comm) {}

            virtual ~DMDABase() {}

            void unit_cube(const SizeType nx, const SizeType ny, const SizeType nz) {
                InputParameters params;
                params.set("nx", nx);
                params.set("ny", ny);
                params.set("nz", nz);

                read(params);
            }

            void read(Input &in) override {
                IntArray dims;

                in.get("nx", dims[0]);
                in.get("ny", dims[1]);
                in.get("nz", dims[2]);

                int n = box_min_.size();
                for (int i = 0; i < n; ++i) {
                    box_min_[i] = 0.0;
                    box_max_[i] = 1.0;
                }

                in.get("min_x", box_min_[0]);
                in.get("min_y", box_min_[1]);
                in.get("min_z", box_min_[2]);

                in.get("max_x", box_max_[0]);
                in.get("max_y", box_max_[1]);
                in.get("max_z", box_max_[2]);

                SizeType n_components = 1;
                in.get("n_components", n_components);

                // FIXME
                create_uniform(comm().get(),
                               dims,
                               box_min_,
                               box_max_,
                               DMDA_ELEMENT_Q1,
                               n_components,
                               raw_type(),
                               DMDA_STENCIL_BOX);
            }

            void set_field_names(const std::vector<std::string> &names) {
                const std::size_t n_names = names.size();
                std::vector<const char *> names_copy(n_names + 1);

                for (std::size_t i = 0; i < n_names; ++i) {
                    names_copy[i] = names[i].c_str();
                }

                names_copy[n_names] = nullptr;

                DMDASetFieldNames(raw_type(), &names_copy[0]);
            }

            void set_field_name(const SizeType &nf, const std::string &name) {
                DMDASetFieldName(raw_type(), nf, name.c_str());
            }

            static void get_corners(DM dm, IntArray &start, IntArray &extent) {
                PetscErrorCode ierr;
                ierr = DMDAGetCorners(dm, &start[0], &start[1], &start[2], &extent[0], &extent[1], &extent[2]);
                assert(ierr == 0);
                UTOPIA_UNUSED(ierr);
            }

            static void get_ghost_corners(DM dm, IntArray &start, IntArray &extent) {
                PetscErrorCode ierr;
                ierr = DMDAGetGhostCorners(dm, &start[0], &start[1], &start[2], &extent[0], &extent[1], &extent[2]);
                assert(ierr == 0);
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

            static void set_dof(DM dm, SizeType n_components) {
                PetscErrorCode ierr;
                ierr = DMDASetDof(dm, n_components);
                assert(ierr == 0);
                UTOPIA_UNUSED(ierr);
            }

            void set_dof(SizeType n_components) { set_dof(raw_type(), n_components); }

            inline void get_corners(IntArray &start, IntArray &extent) const { get_corners(raw_type(), start, extent); }

            inline void get_ghost_corners(IntArray &start, IntArray &extent) const {
                get_ghost_corners(raw_type(), start, extent);
            }

            inline SizeType get_dof() const { return get_dof(raw_type()); }

            static void uniform_refine(DM dm, DM *fine_dm) {
                MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(dm));
                DMRefine(dm, comm, fine_dm);

                // This does not transfer automatically for some reason
                DMDAElementType elem_type;
                DMDAGetElementType(dm, &elem_type);
                DMDASetElementType(dm, elem_type);
            }

            std::unique_ptr<DMDABase> uniform_refine() {
                std::unique_ptr<DMDABase> fine = utopia::make_unique<DMDABase>(this->comm());
                uniform_refine(raw_type(), &fine->raw_type());

                fine->box_min_.copy(box_min_);
                fine->box_max_.copy(box_max_);
                return fine;
            }

            std::unique_ptr<DMDABase> clone() {
                std::unique_ptr<DMDABase> cloned = utopia::make_unique<DMDABase>(this->comm());
                DMClone(raw_type(), &cloned->raw_type());
                cloned->box_min_.copy(box_min_);
                cloned->box_max_.copy(box_max_);
                return cloned;
            }

            inline DMDAElementType elem_type() const {
                DMDAElementType ret;
                DMDAGetElementType(raw_type(), &ret);
                return ret;
            }

            static void get_dims(DM dm, IntArray &dims) {
                PetscErrorCode ierr;

                ierr = DMDAGetInfo(dm,
                                   nullptr,
                                   &dims[0],
                                   &dims[1],
                                   &dims[2],
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
                UTOPIA_UNUSED(ierr);
            }

            inline void get_dims(IntArray &dims) const { get_dims(raw_type(), dims); }

            static void create(const MPI_Comm &comm,
                               const IntArray &dims,
                               const SizeType n_components,
                               DMDAStencilType stencil_type,
                               DM &dm) {
                SizeType dim = 0;

                for (auto d : dims) {
                    dim += d != 0;
                }

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

                SizeType dim = 0;

                for (auto d : dims) {
                    dim += d != 0;
                }

                if (dim > 1) {
                    min_y = box_min[1];
                    max_y = box_max[1];
                }

                if (dim > 2) {
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

            inline SizeType elements_x_cell() const {
                int exc = 1;
                if (elem_type() == DMDA_ELEMENT_P1) {
                    int dim = this->get_dimension();
                    if (dim == 2) {
                        exc = 2;
                    } else if (dim == 3) {
                        exc = 6;
                    } else {
                        assert(false);
                    }
                }

                return exc;
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

            std::shared_ptr<Elements> elements_;
            Point box_min_;
            Point box_max_;

            inline const Point &box_min() const { return box_min_; }
            inline const Point &box_max() const { return box_max_; }
        };

        StructuredGrid::StructuredGrid(const Communicator &comm) : impl_(utopia::make_unique<DMDABase>(comm)) {}
        StructuredGrid::~StructuredGrid() = default;

        void StructuredGrid::unit_cube(const SizeType nx, const SizeType ny, const SizeType nz) {
            // Add one to be conforming to the other backends
            impl_->unit_cube(nx + 1, ny + 1, nz + 1);
        }

        StructuredGrid::Communicator &StructuredGrid::comm() { return impl_->comm(); }
        const StructuredGrid::Communicator &StructuredGrid::comm() const { return impl_->comm(); }

        std::unique_ptr<StructuredGrid> StructuredGrid::uniform_refine() const {
            auto fine = utopia::make_unique<StructuredGrid>(comm());
            fine->impl_ = std::move(impl_->uniform_refine());
            return fine;
        }

        std::unique_ptr<StructuredGrid> StructuredGrid::clone(const SizeType &n_components) const {
            auto cloned = utopia::make_unique<StructuredGrid>(comm());
            cloned->impl_ = std::move(impl_->clone());
            // cloned->impl_->set_dof(n_components);
            return cloned;
        }

        std::unique_ptr<StructuredGrid> StructuredGrid::clone() const { return clone(this->impl_->get_dof()); }

        void StructuredGrid::read(Input &in) { impl_->read(in); }

        void StructuredGrid::set_field_name(const SizeType &nf, const std::string &name) {
            impl_->set_field_name(nf, name);
        }

        void StructuredGrid::set_field_names(const std::vector<std::string> &names) { impl_->set_field_names(names); }

        // void StructuredGrid::nodes(const SizeType &idx, NodeIndex &nodes) const {
        //     assert(elements_ && "init_elements() has to be called before any access to element data");
        //     nodes = impl_->elements_->nodes(idx);
        // }

        // void StructuredGrid::nodes_local(const SizeType &idx, NodeIndex &nodes) const {
        //     assert(elements_ && "init_elements() has to be called before any access to element data");
        //     nodes = impl_->elements_->nodes_local(idx);
        // }

        // bool StructuredGrid::on_boundary(const SizeType &elem_idx) const {
        //     NodeIndex idx;
        //     nodes(elem_idx, idx);

        //     for (auto i : idx) {
        //         if (this->impl_->is_node_on_boundary(i)) {
        //             return true;
        //         }
        //     }

        //     return false;
        // }

        StructuredGrid::View StructuredGrid::view() const {
            SizeType dof_range_begin, dof_range_end;
            IntArray dims, corners_begin, corners_extent, ghost_corners_begin, ghost_corners_extent;

            impl_->dof_ownership_range(dof_range_begin, dof_range_end);
            const SizeType dim = impl_->get_dimension();
            impl_->get_dims(dims);
            impl_->get_corners(corners_begin, corners_extent);
            impl_->get_ghost_corners(ghost_corners_begin, ghost_corners_extent);
            SizeType n_components = impl_->get_dof();
            auto &box_min = impl_->box_min();
            auto &box_max = impl_->box_max();

            std::vector<Scalar> v_box_min(dim), v_box_max(dim);
            std::vector<SizeType> v_dims(dim), v_corners_begin(dim), v_corners_extent(dim), v_ghost_corners_begin(dim),
                v_ghost_corners_extent(dim);

            for (SizeType d = 0; d < dim; ++d) {
                v_dims[d] = dims[d];
                v_corners_begin[d] = corners_begin[d];
                v_corners_extent[d] = corners_extent[d];
                v_ghost_corners_begin[d] = ghost_corners_begin[d];
                v_ghost_corners_extent[d] = ghost_corners_extent[d];
                v_box_min[d] = box_min[d];
                v_box_max[d] = box_max[d];
            }

            View v(v_dims,
                   v_corners_begin,
                   v_corners_extent,
                   v_ghost_corners_begin,
                   v_ghost_corners_extent,
                   v_box_min,
                   v_box_max,
                   n_components,
                   impl_->elements_x_cell());

            v.set_dof_range_begin(dof_range_begin);
            v.set_dof_range_end(dof_range_end);
            return v;
        }

        void StructuredGrid::describe(std::ostream &os) const {
            auto v = view();
            os << "n_elements      : " << v.n_elements() << '\n';
            os << "n_nodes         : " << v.n_nodes() << '\n';
            os << "dim             : " << v.dim() << '\n';
            os << "elements_x_cell : " << v.elements_x_cell() << '\n';
        }

        StructuredGrid::SizeType StructuredGrid::n_nodes() const {
            IntArray dims;
            impl_->get_dims(dims);

            SizeType ret = 1;
            for (auto d : dims) {
                if (d > 0) ret *= d;
            }
            return ret;
        }

    }  // namespace petsc
}  // namespace utopia
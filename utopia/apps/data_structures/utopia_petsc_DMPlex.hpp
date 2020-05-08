#ifndef UTOPIA_PETSC_DM_PLEX_HPP
#define UTOPIA_PETSC_DM_PLEX_HPP

#include "utopia_StructuredGrid.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_petsc_DM.hpp"

#include "utopia_Algorithms.hpp"

#include <petscdmplex.h>
#include <cassert>

namespace utopia {

    template <class Point, class IntArray>
    class PetscDMPlex;

    template <class Point, class IntArray>
    class PetscDMPlex : public PetscDMBase {
    public:
        // using Super = utopia::StructuredGrid<Point, IntArray>;

        using SizeType = typename Traits<IntArray>::ValueType;
        using Scalar = typename Traits<Point>::Scalar;
        using NodeIndex = utopia::ArrayView<SizeType>;

        // class Elements {
        // public:
        //     Elements(DM dm)
        //     : dm_(dm)
        //     {
        //         DMPlexGetElements(dm_, &n_local_elem_, &n_nodes_x_elem_, &local_elem_);
        //     }

        //     ~Elements()
        //     {
        //         DMPlexRestoreElements(dm_, &n_local_elem_, &n_nodes_x_elem_, &local_elem_);
        //     }

        //     inline SizeType local_size() const
        //     {
        //         return n_local_elem_;
        //     }

        //     inline NodeIndex nodes_local(const SizeType &local_elem_idx) const
        //     {
        //         return NodeIndex(
        //             &local_elem_[local_elem_idx * n_nodes_x_elem_],
        //             n_nodes_x_elem_
        //         );
        //     }

        // private:
        //     DM dm_;
        //     SizeType n_local_elem_, n_nodes_x_elem_;
        //     const SizeType *local_elem_;
        // };

        // std::unique_ptr<Elements> make_elements() const
        // {
        //     return utopia::make_unique<Elements>(this->raw_type());
        // }

        void wrap(DM &dm, const bool delegate_ownership) override {
            PetscDMBase::wrap(dm, delegate_ownership);
            init_from_dm(dm);
        }

        PetscDMPlex(const PetscCommunicator &comm) : PetscDMBase(comm) {}

        PetscDMPlex(DM &dm, const bool delegate_ownership) { wrap(dm, delegate_ownership); }

        bool read(const Path &path, const bool interpolate = false) {
            PetscErrorCode ierr = 0;

            const auto ext = path.extension();

            if (ext == "e") {
                this->destroy_dm();

                PetscBool p_inter = interpolate ? PETSC_TRUE : PETSC_FALSE;
                ierr = DMPlexCreateExodusFromFile(comm().get(), path.c_str(), p_inter, &this->raw_type());
                assert(ierr == 0);
            } else {
                return false;
            }

            return ierr == 0;
        }

        void read(Input &in) override {
            this->destroy_dm();

            init_default();

            PetscErrorCode ierr = 0;
            std::string type = "box";
            SizeType dim = 2;
            bool simplex = false;
            bool interpolate = false;
            Scalar lower[3] = {0, 0, 0};
            Scalar upper[3] = {1, 1, 1};
            SizeType faces[3] = {2, 2, 2};
            Scalar refinement_limit = 0.0;

            in.get("type", type);
            in.get("dim", dim);
            in.get("simplex", simplex);
            in.get("interpolate", interpolate);

            in.get("x_min", lower[0]);
            in.get("y_min", lower[1]);
            in.get("z_min", lower[2]);

            in.get("x_max", upper[0]);
            in.get("y_max", upper[1]);
            in.get("z_max", upper[2]);

            in.get("nx", faces[0]);
            in.get("ny", faces[1]);
            in.get("nz", faces[2]);
            in.get("refinement_limit", refinement_limit);

            if (type == "box") {
                ierr = DMPlexCreateBoxMesh(this->comm().get(),
                                           dim,
                                           simplex ? PETSC_TRUE : PETSC_FALSE,
                                           faces,
                                           lower,
                                           upper,
                                           nullptr,
                                           interpolate ? PETSC_TRUE : PETSC_FALSE,
                                           &this->raw_type());
                assert(ierr == 0);
            } else if (type == "file") {
                std::string path;
                in.get("path", path);

                if (path.empty()) {
                    std::cerr << "[Error] empty path" << std::endl;
                    return;
                }

                if (!this->read(path.c_str())) {
                    std::cerr << "[Error] bad path: " << path << std::endl;
                    return;
                }

                if (interpolate) {
                    DMPlexInterpolate(this->raw_type(), &this->raw_type());

                    // DM dm = nullptr;
                    // DMPlexInterpolate(this->raw_type(), &dm);
                    // DMPlexCopyCoordinates(this->raw_type(), dm);
                    // this->destroy_dm();
                    // this->raw_type() = dm;
                }
            }

            ierr = PetscObjectSetName((PetscObject)raw_type(), "Mesh");

            partition(refinement_limit);
        }

        // const SizeType n = this->dims().size();
        // assert(n > 0); //IMPLEMENT ME for dynamically sized arrays

        // const char coord_names[3] = {'x', 'y', 'z'};

        // std::string n_str   = "n?";
        // std::string min_str = "?_min";
        // std::string max_str = "?_max";

        // for(SizeType d = 0; d < n; ++d) {
        //     const char coord = coord_names[d];

        //     n_str[1]   = coord;
        //     min_str[0] = coord;
        //     max_str[0] = coord;

        //     in.get(n_str,   this->dims()[d]);
        //     in.get(min_str, this->box_min()[d]);
        //     in.get(max_str, this->box_max()[d]);
        // }

        // this->destroy_dm();
        // create_uniform(
        //     comm().get(),
        //     this->dims(),
        //     this->box_min(),
        //     this->box_max(),
        //     type_override_,
        //     this->n_components(),
        //     this->raw_type()
        // );

        /// re-initialize mirror just to be safe (it is cheap anyway)
        //     update_mirror();
        // }

        void update_mirror() { init_from_dm(raw_type()); }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::unique_ptr<PetscDMPlex> uniform_refine() const {
            auto fine = utopia::make_unique<PetscDMPlex>(comm());
            PetscDMBase::refine(raw_type(), comm().get(), fine->raw_type());

            // This does not transfer automatically for some reason
            // DMPlexElementType elem_type;
            // DMPlexGetElementType(raw_type(), &elem_type);
            // DMPlexSetElementType(fine->raw_type(), elem_type);
            fine->update_mirror();

            return std::move(fine);
        }

        // std::unique_ptr<PetscDMPlex> clone(const SizeType &n_components) const
        // {
        //     auto cloned = utopia::make_unique<PetscDMPlex>();
        //     cloned->copy(*this);
        //     cloned->set_n_components(n_components);
        //     cloned->type_override_ = type_override_;
        //     cloned->init_from_mirror();
        //     return std::move(cloned);
        // }

        inline SizeType dim() const { return PetscDMBase::get_dimension(this->raw_type()); }

    private:
        // DMPlexElementType type_override_;

        void init_default() {
            // device::fill(10, this->dims());
            // device::fill(0,  this->box_min());
            // device::fill(1,  this->box_max());
            // this->set_n_components(1);
        }

        void init_from_mirror() {
            // this->destroy_dm();
            // create_uniform(
            //     comm().get(),
            //     this->dims(),
            //     this->box_min(),
            //     this->box_max(),
            //     type_override_,
            //     this->n_components(),
            //     this->raw_type()
            // );
        }

        void init_from_dm(DM dm) {
            MPI_Comm mpi_comm = PetscObjectComm((PetscObject)dm);
            comm().set(mpi_comm);

            PetscInt dof_range_begin, dof_range_end;
            PetscDMBase::dof_ownership_range(dm, dof_range_begin, dof_range_end);

            const SizeType dim = PetscDMBase::get_dimension(dm);
            // assert(dim == this->dim());

            // this->set_dof_range_begin(dof_range_begin);
            // this->set_dof_range_end(dof_range_end);

            // get_dims(dm, this->dims());
            // get_corners(dm, this->corners_begin(), this->corners_extent());
            // get_ghost_corners(dm, this->ghost_corners_begin(), this->ghost_corners_extent());
            // this->set_n_components( get_dof(dm) );

            // DMPlexElementType elem_type;
            // DMPlexGetElementType(dm, &elem_type);

            // if(elem_type == DMPlex_ELEMENT_P1) {
            //     if(this->dim() == 2) {
            //         this->set_elements_x_cell(2);
            //     } else if(this->dim() == 3) {
            //         this->set_elements_x_cell(6);
            //     } else {
            //         assert(false);
            //     }
            // }
        }

        void partition(const Scalar &refinement_limit) {
            PetscErrorCode ierr = 0;
            PetscPartitioner part;
            DM refined_mesh = NULL;
            DM distributed_mesh = NULL;

            auto *dm = &raw_type();

            /* Refine mesh using a volume constraint */
            if (refinement_limit > 0.0) {
                ierr = DMPlexSetRefinementLimit(*dm, refinement_limit);
                assert(ierr == 0);
                ierr = DMRefine(*dm, this->comm().get(), &refined_mesh);
                assert(ierr == 0);
                if (refined_mesh) {
                    const char *name;

                    ierr = PetscObjectGetName((PetscObject)*dm, &name);
                    assert(ierr == 0);
                    ierr = PetscObjectSetName((PetscObject)refined_mesh, name);
                    assert(ierr == 0);
                    ierr = DMDestroy(dm);
                    assert(ierr == 0);
                    *dm = refined_mesh;
                }
            }
            /* Distribute mesh over processes */

            ierr = DMPlexGetPartitioner(*dm, &part);
            assert(ierr == 0);
            ierr = PetscPartitionerSetFromOptions(part);
            assert(ierr == 0);
            ierr = DMPlexDistribute(*dm, 0, NULL, &distributed_mesh);
            assert(ierr == 0);
            if (distributed_mesh) {
                ierr = DMDestroy(dm);
                assert(ierr == 0);
                *dm = distributed_mesh;
            }
        }
    };  // namespace utopia

}  // namespace utopia

#endif  // UTOPIA_PETSC_DM_PLEX_HPP

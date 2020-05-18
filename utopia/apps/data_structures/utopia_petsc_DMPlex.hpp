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

    template <class Point_, class IntArray>
    class PetscDMPlex : public PetscDMBase {
    public:
        // using Super = utopia::StructuredGrid<Point, IntArray>;

        using Point = Point_;
        using SizeType = typename Traits<IntArray>::ValueType;
        using Scalar = typename Traits<Point>::Scalar;
        using NodeIndex = utopia::ArrayView<SizeType>;

        using Device = utopia::Device<PETSC>;
        using Comm = utopia::PetscCommunicator;

        void wrap(DM &dm, const bool delegate_ownership) override {
            PetscDMBase::wrap(dm, delegate_ownership);
            init_from_dm(dm);
        }

        PetscDMPlex(const PetscCommunicator &comm) : PetscDMBase(comm), interpolated_(false), simplex_(false) {}

        // FIXME interpolated_
        PetscDMPlex(DM &dm, const bool delegate_ownership) : interpolated_(false), simplex_(false) {
            wrap(dm, delegate_ownership);
        }

        bool read(const Path &path, const bool interpolate = false) {
            PetscErrorCode ierr = 0;

            const auto ext = path.extension();
            PetscBool p_inter = interpolate ? PETSC_TRUE : PETSC_FALSE;

            // if (ext == "e") {
            //     this->destroy_dm();
            //     ierr = DMPlexCreateExodusFromFile(comm().get(), path.c_str(), p_inter, &this->raw_type());
            //     assert(ierr == 0);
            // } else {
            ierr = DMPlexCreateFromFile(comm().get(), path.c_str(), p_inter, &this->raw_type());
            assert(ierr == 0);
            // }

            return ierr == 0;
        }

        void read(Input &in) override {
            this->destroy_dm();

            init_default();

            PetscErrorCode ierr = 0;
            std::string type = "box";
            SizeType dim = 2;

            Scalar lower[3] = {0, 0, 0};
            Scalar upper[3] = {1, 1, 1};
            SizeType faces[3] = {2, 2, 2};
            Scalar refinement_limit = 0.0;

            in.get("type", type);
            in.get("dim", dim);
            in.get("simplex", simplex_);
            in.get("interpolate", interpolated_);

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
                                           simplex_ ? PETSC_TRUE : PETSC_FALSE,
                                           faces,
                                           lower,
                                           upper,
                                           nullptr,
                                           interpolated_ ? PETSC_TRUE : PETSC_FALSE,
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

                if (interpolated_) {
                    DMPlexInterpolate(this->raw_type(), &this->raw_type());

                    // DM dm = nullptr;
                    // DMPlexInterpolate(this->raw_type(), &dm);
                    // DMPlexCopyCoordinates(this->raw_type(), dm);
                    // this->destroy_dm();
                    // this->raw_type() = dm;
                }
            }

            bool extrude = false;
            in.get("extrude", extrude);

            if (extrude) {
                SizeType extrude_layers = 1;
                Scalar extrude_height = 1;
                bool extrude_sort = false;

                in.get("extrude_layers", extrude_layers);
                in.get("extrude_height", extrude_height);
                in.get("extrude_sort", extrude_sort);

                DM extrude_dm = nullptr;
                DMPlexExtrude(raw_type(),
                              extrude_layers,
                              extrude_height,
                              extrude_sort ? PETSC_TRUE : PETSC_FALSE,
                              interpolated_ ? PETSC_TRUE : PETSC_FALSE,
                              &extrude_dm);

                this->destroy_dm();
                this->raw_type() = extrude_dm;
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

        inline SizeType n_components() const {
            PetscSection section = nullptr;
            DMGetLocalSection(this->raw_type(), &section);
            SizeType num_comp;
            PetscSectionGetFieldComponents(section, 0, &num_comp);
            return num_comp;
        }

        inline SizeType component(const SizeType &idx) const {
            const SizeType nc = n_components();
            return nc == 1 ? 0 : idx % nc;
        }

        bool is_node_on_boundary(const SizeType &idx, const SideSet::BoundaryIdType b_id) const {
            // DMLabelGetValue(DMLabel label, PetscInt point, PetscInt * value)
            assert(false);
            return false;
        }

        inline SizeType dim() const { return PetscDMBase::get_dimension(this->raw_type()); }

        inline SizeType n_local_elements() const {
            SizeType ret;
            DMPlexGetHeightStratum(this->raw_type(), 0, nullptr, &ret);
            return ret;
        }

        inline Range element_range() const { return Range(0, n_local_elements()); }

        inline bool interpolated() const {
            // DMPlexInterpolatedFlag ret;
            // DMPlexIsInterpolated(this->raw_type(), &ret);
            // return DMPLEX_INTERPOLATED_INVALID != ret;
            return interpolated_;
        }

        inline void set_num_fields(const SizeType &num) { DMSetNumFields(this->raw_type(), num); }

        void create_section(const SizeType *num_comp, const SizeType *num_dofs) {
            PetscSection section = nullptr;
            DMPlexCreateSection(
                this->raw_type(), nullptr, num_comp, num_dofs, 0, nullptr, nullptr, nullptr, nullptr, &section);

            DMSetLocalSection(this->raw_type(), section);
            PetscSectionDestroy(&section);
        }

        void set_field_name(const SizeType &field, const std::string &name) {
            PetscSection section = nullptr;
            DMGetLocalSection(this->raw_type(), &section);
            PetscSectionSetFieldName(section, field, name.c_str());
        }

        inline void set_up() {
            DMSetUp(this->raw_type());
            init_coords();
        }

        inline void transform(const SizeType &cell, const Point &ref, Point &physical) const {
            DMPlexReferenceToCoordinates(this->raw_type(), cell, 1, &ref[0], &physical[0]);
        }

        template <typename IntArrayT>
        inline void nodes(const SizeType &cell, IntArrayT &nodes) const {
            SizeType num_points = 0;
            SizeType *points = nullptr;

            DMPlexGetTransitiveClosure(this->raw_type(), cell, PETSC_TRUE, &num_points, &points);

            for (SizeType i = 1; i < num_points; ++i) {
                nodes[i - 1] = points[i * 2];
            }

            DMPlexRestoreTransitiveClosure(this->raw_type(), cell, PETSC_TRUE, &num_points, &points);
        }

        template <typename IntArrayT>
        inline void nodes_local(const SizeType &cell, IntArrayT &nodes) const {
            SizeType num_points = 0;
            SizeType *points = nullptr;

            DMPlexGetTransitiveClosure(this->raw_type(), cell, PETSC_TRUE, &num_points, &points);

            SizeType end_dummy;
            for (SizeType i = 1; i < num_points; ++i) {
                DMPlexGetPointLocal(this->raw_type(), points[i * 2], &nodes[i - 1], &end_dummy);
            }

            DMPlexRestoreTransitiveClosure(this->raw_type(), cell, PETSC_TRUE, &num_points, &points);
        }

        template <typename IntArrayT>
        inline void fields_local(const SizeType &cell, const SizeType &field, IntArrayT &nodes) {
            SizeType num_points = 0;
            SizeType *points = nullptr;

            DMPlexGetTransitiveClosure(this->raw_type(), cell, PETSC_TRUE, &num_points, &points);

            SizeType end_dummy;
            for (SizeType i = 1; i < num_points; ++i) {
                DMPlexGetPointLocalField(this->raw_type(), points[i * 2], field, &nodes[i - 1], &end_dummy);
            }

            DMPlexRestoreTransitiveClosure(this->raw_type(), cell, PETSC_TRUE, &num_points, &points);
        }

        template <typename IntArrayT>
        inline void cone(const SizeType &cell, IntArrayT &nodes) {
            const SizeType *cone = nullptr;
            DMPlexGetCone(this->raw_type(), cell, &cone);

            SizeType num_faces = 0;
            DMPlexGetConeSize(this->raw_type(), cell, &num_faces);

            for (SizeType i = 0; i < num_faces; ++i) {
                nodes[i] = cone[i];
            }
        }

        inline void simplex(const bool simplex) { simplex_ = simplex; }

        void describe(std::ostream &os = std::cout) const {
            SizeType nl, num_cs, num_vs, num_fs, num_marker, num_bd;

            DMGetNumLabels(this->raw_type(), &nl);

            os << "nl:  " << nl << std::endl;

            for (SizeType i = 0; i < nl; ++i) {
                const char *label_name;
                DMGetLabelName(this->raw_type(), i, &label_name);

                os << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
                os << i << ") " << label_name << std::endl;

                // DMGetLabelSize(this->raw_type(), "Cell Sets", &num_cs);
                // DMGetLabelSize(this->raw_type(), "Vertex Sets", &num_vs);
                // DMGetLabelSize(this->raw_type(), "Face Sets", &num_fs);
                // DMGetLabelSize(this->raw_type(), "marker", &num_marker);

                // os << "num_cs: " << num_cs << " num_vs: " << num_vs << " num_fs: " << num_fs
                //    << " num_marker: " << num_marker << std::endl;

                // DMGetNumBoundary(this->raw_type(), &num_bd);

                // os << "num_bd: " << num_bd << std::endl;

                DMLabel label;
                // DMGetLabel(this->raw_type(), "marker", &label);
                // DMGetLabel(this->raw_type(), "Cell Sets", &label);
                DMGetLabel(this->raw_type(), label_name, &label);

                if (!label) return;

                for (SizeType s = 0; s < 4; ++s) {
                    os << "---------------------------\n";
                    os << "Stratum " << s << std::endl;

                    SizeType n;
                    DMLabelGetStratumSize(label, s, &n);

                    os << "DMLabelGetStratumSize " << n << std::endl;

                    if (n == 0) continue;

                    DMLabelGetNumValues(label, &n);

                    os << "DMLabelGetNumValues " << n << std::endl;

                    SizeType start, end;
                    DMPlexGetHeightStratum(this->raw_type(), s, &start, &end);

                    os << "[" << start << ", " << end << ")" << std::endl;

                    for (SizeType i = start; i < end; ++i) {
                        SizeType value;
                        DMLabelGetValue(label, i, &value);
                        os << value << " ";
                    }

                    os << std::endl;
                }

                os << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
            }
        }

    private:
        // DMPlexElementType type_override_;
        bool interpolated_;
        bool simplex_;
        PetscVector coords_;

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

        void init_coords() {
            Vec coords;
            // DMGetCoordinates(dmplex.raw_type(), &coords.raw_type());
            DMGetCoordinatesLocal(this->raw_type(), &coords);
            // We do not need to delete the memory
            coords_.wrap(coords);
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

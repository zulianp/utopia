#ifndef UTOPIA_PETSC_DMDA_HPP
#define UTOPIA_PETSC_DMDA_HPP

#include "utopia_petsc.hpp"
#include "utopia_petsc_DMForwardDeclarations.hpp"

#include "utopia_Algorithms.hpp"
#include "utopia_Traits.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_AABB.hpp"
#include "utopia_mesh_StructuredGrid.hpp"

#include <cassert>

namespace utopia {

    template <>
    class Traits<petsc::StructuredGrid> : public Traits<PetscVector> {
    public:
        using SideSet = utopia::mesh::SideSet;
        using AABB = utopia::AABB<std::vector<Traits<PetscVector>::Scalar>>;
    };

    namespace petsc {

        class DMDABase;

        class StructuredGrid : public Configurable, public Describable, public Traits<StructuredGrid> {
        public:
            using Traits = utopia::Traits<StructuredGrid>;
            using Communicator = Traits::Communicator;
            using SizeType = Traits::SizeType;
            using Scalar = Traits::Scalar;
            using NodeIndex = utopia::ArrayView<const SizeType>;
            using IntArray = utopia::ArrayView<SizeType, 3>;
            using Point = utopia::StaticVector<Scalar, 3>;
            using View = mesh::StructuredGrid<std::vector<Scalar>, std::vector<SizeType>>;
            using ElemToNodeIndex = utopia::ArrayView<const SizeType, DYNAMIC_SIZE, DYNAMIC_SIZE>;

            // bool read(const Path &path);
            bool write(const Path &path) const;

            // FIXME make this work also without static-sizes
            // constexpr static typename SideSets::Sides sides() { return SideSets::sides(); }

            // void nodes(const SizeType &idx, NodeIndex &nodes) const;
            // void nodes_local(const SizeType &idx, NodeIndex &nodes) const;
            // bool on_boundary(const SizeType &elem_idx) const;

            void set_field_name(const SizeType &nf, const std::string &name);
            void set_field_names(const std::vector<std::string> &names);

            void read(Input &in) override;
            void unit_cube(const SizeType nx, const SizeType ny, const SizeType nz);

            void update_mirror();

            void bounding_box(AABB &output) const;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            std::unique_ptr<StructuredGrid> uniform_refine() const;
            std::unique_ptr<StructuredGrid> clone(const SizeType &n_components) const;
            std::unique_ptr<StructuredGrid> clone() const;

            void describe(std::ostream &os) const override;

            Communicator &comm();
            const Communicator &comm() const;

            View view() const;

            StructuredGrid(const Communicator &comm);
            ~StructuredGrid();

            SizeType n_elements() const;
            SizeType n_nodes() const;

            SizeType n_local_elements() const;
            SizeType n_local_nodes() const;

            DMDABase &dm() const;

            ElemToNodeIndex elem_to_node_index() const;
            ElemToNodeIndex elem_to_local_node_index() const;

            int spatial_dimension() const;

            void box(const AABB &box,
                     const SizeType &nx,
                     const SizeType &ny,
                     const SizeType &nz,
                     const std::string &elem_type);

        private:
            std::unique_ptr<DMDABase> impl_;
        };

    }  // namespace petsc
}  // namespace utopia

#endif  // UTOPIA_PETSC_DMDA_HPP

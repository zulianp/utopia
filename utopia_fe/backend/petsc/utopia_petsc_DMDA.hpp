#ifndef UTOPIA_PETSC_DMDA_HPP
#define UTOPIA_PETSC_DMDA_HPP

#include "utopia_petsc.hpp"
#include "utopia_petsc_DMForwardDeclarations.hpp"

#include "utopia_Algorithms.hpp"
#include "utopia_Traits.hpp"
// #include "utopia_StructuredGrid.hpp"
#include "utopia_make_unique.hpp"

// #include "utopia_petsc_FE.hpp"

#include <cassert>

namespace utopia {

    template <>
    class Traits<petsc::StructuredGrid> : public Traits<PetscVector> {};

    namespace petsc {

        class StructuredGrid : public Configurable, public Describable {
        public:
            using Traits = utopia::Traits<StructuredGrid>;
            using Communicator = Traits::Communicator;
            using SizeType = Traits::SizeType;
            using Scalar = Traits::Scalar;
            using NodeIndex = utopia::ArrayView<const SizeType>;
            using IntArray = utopia::ArrayView<SizeType, 3>;
            using Point = utopia::StaticVector<Scalar, 3>;

            // FIXME make this work also without static-sizes
            // constexpr static typename SideSets::Sides sides() { return SideSets::sides(); }

            void nodes(const SizeType &idx, NodeIndex &nodes) const;
            void nodes_local(const SizeType &idx, NodeIndex &nodes) const;
            bool on_boundary(const SizeType &elem_idx) const;

            void set_field_name(const SizeType &nf, const std::string &name);
            void set_field_names(const std::vector<std::string> &names);

            void read(Input &in) override;

            void update_mirror();
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            std::unique_ptr<StructuredGrid> uniform_refine() const;
            std::unique_ptr<StructuredGrid> clone(const SizeType &n_components) const;
            std::unique_ptr<StructuredGrid> clone() const;

            void describe(std::ostream &os) const override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace petsc
}  // namespace utopia

#endif  // UTOPIA_PETSC_DMDA_HPP

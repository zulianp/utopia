#ifndef UTOPIA_LIBMESH_LEGACY_HPP
#define UTOPIA_LIBMESH_LEGACY_HPP

// New
#include "utopia_libmesh_FunctionSpace_new.hpp"

// Legacy
#include "utopia_ProductFunctionSpace.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

#include <memory>

namespace utopia {
    namespace libmesh {

        using LegacyFunctionSpace = utopia::LibMeshFunctionSpace;
        using LegacyProductFunctionSpace = utopia::ProductFunctionSpace<LegacyFunctionSpace>;

        inline std::shared_ptr<LegacyProductFunctionSpace> make_legacy(utopia::libmesh::FunctionSpace &space) {
            auto legacy_space = std::make_shared<LegacyProductFunctionSpace>();

            for (int s = 0; s < space.n_subspaces(); ++s) {
                legacy_space->add_subspace(std::make_shared<LegacyFunctionSpace>(
                    make_ref(space.raw_type()), space.system_id(), space[s].subspace_id()));
            }

            return legacy_space;
        }

    }  // namespace libmesh
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_LEGACY_HPP
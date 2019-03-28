#ifndef UTOPIA_LIBMESH_TREE_NAVIGATOR_HPP
#define UTOPIA_LIBMESH_TREE_NAVIGATOR_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_Types.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_FunctionalTraits.hpp"
#include "utopia_libmesh_Utils.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_Traverse.hpp"
#include "utopia_ProductFunctionSpace.hpp"
#include "utopia_libmesh_TreeNavigator.hpp"
#include "utopia_FindSpace.hpp"

namespace utopia {

    template<class Expr>
    static std::shared_ptr<LibMeshFunctionSpace> find_any_space(const Expr &expr)
    {
        std::shared_ptr<LibMeshFunctionSpace> space_ptr = test_space<LibMeshFunctionSpace>(expr);
        if(space_ptr) return space_ptr;

        auto prod_space_ptr = test_space<ProductFunctionSpace<LibMeshFunctionSpace>>(expr);
        if(prod_space_ptr) {
            space_ptr = prod_space_ptr->subspace_ptr(0);
        }

        if(space_ptr) return space_ptr;

        return find_trial_space<LibMeshFunctionSpace>(expr);

        // space_ptr = trial_space<LibMeshFunctionSpace>(expr);
        // if(space_ptr) return space_ptr;

        // prod_space_ptr = trial_space<ProductFunctionSpace<LibMeshFunctionSpace>>(expr);

        // if(prod_space_ptr) {
        // 	space_ptr = prod_space_ptr->subspace_ptr(0);
        // }

        // return space_ptr;
    }
}

#endif //UTOPIA_LIBMESH_TREE_NAVIGATOR_HPP

#include "utopia_libmesh_Deprecated.hpp"

#include "moonolith_vector.hpp"
#include "utopia_intersector.hpp"
#include "utopia_libmesh_Transform.hpp"
#include "utopia_libmesh_Utils.hpp"

#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {
    //

    // template <int Dim>
    // void ConvertFunctionSpace<LibMeshFunctionSpace, moonolith::FunctionSpace<moonolith::Mesh<double, Dim>>>::apply(
    //     const LibMeshFunctionSpace &in_space,
    //     moonolith::FunctionSpace<moonolith::Mesh<double, Dim>> &out) {
    //     int var_num = in_space.subspace_id();
    //     auto &in = in_space.mesh();
    //     auto &dof_map = in_space.dof_map();
    //     apply(in, dof_map, var_num, out);
    // }

    // template class ConvertFunctionSpace<LibMeshFunctionSpace, moonolith::FunctionSpace<moonolith::Mesh<double, 1>>>;
    // template class ConvertFunctionSpace<LibMeshFunctionSpace, moonolith::FunctionSpace<moonolith::Mesh<double, 2>>>;
    // template class ConvertFunctionSpace<LibMeshFunctionSpace, moonolith::FunctionSpace<moonolith::Mesh<double, 3>>>;
}  // namespace utopia

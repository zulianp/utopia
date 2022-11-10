#ifndef UTOPIA_LIBMESH_DEPRECATED_HPP
#define UTOPIA_LIBMESH_DEPRECATED_HPP

#include "moonolith_function_space.hpp"
#include "utopia_libmesh_Transform.hpp"

#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_moonolith_libmesh_Convert.hpp"

#include "utopia_intersector.hpp"
#include "utopia_libmesh_QMortar.hpp"

namespace utopia {

    // void compute_side_normal(const int dim, const libMesh::Elem &side, libMesh::Point &n);

    // int order_for_l2_integral(const int dim,
    //                           const libMesh::Elem &master_el,
    //                           const int master_order,
    //                           const libMesh::Elem &slave_el,
    //                           const int slave_order);

    // void print(const libMesh::QBase &ir, std::ostream &os = std::cout);
    // double sum_of_weights(const libMesh::QBase &ir);
    // double sum(const libMesh::DenseMatrix<libMesh::Real> &mat);

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_DEPRECATED_HPP

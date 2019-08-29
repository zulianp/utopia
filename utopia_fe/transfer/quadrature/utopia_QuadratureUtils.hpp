#ifndef UTOPIA_QUADRATURE_UTILS_HPP
#define UTOPIA_QUADRATURE_UTILS_HPP

#include <libmesh/quadrature_gauss.h>
#include <libmesh/elem.h>
#include <libmesh/fe_type.h>

#include <memory>

namespace utopia {

    class QuadratureUtils {
    public:
        static std::shared_ptr<libMesh::QBase> nodal_quad_points(
            const libMesh::Elem &elem);
    };
}

#endif //UTOPIA_QUADRATURE_UTILS_HPP

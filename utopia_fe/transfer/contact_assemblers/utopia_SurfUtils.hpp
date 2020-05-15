#ifndef UTOPIA_SURF_UTILS_HPP
#define UTOPIA_SURF_UTILS_HPP

#include "libmesh/dense_matrix.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"

namespace utopia {
    class SurfUtils {
    public:
        using Elem = libMesh::Elem;
        using FEType = libMesh::FEType;
        using Real = libMesh::Real;
        using Point = libMesh::Point;

        static void avg_normal(const Elem &trial, const FEType &type, Point &normal);
    };
}  // namespace utopia

#endif  // UTOPIA_SURF_UTILS_HPP

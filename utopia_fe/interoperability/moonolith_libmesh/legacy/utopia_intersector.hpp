#ifndef UTOPIA_INTERSECTOR_HPP
#define UTOPIA_INTERSECTOR_HPP

#include <libmesh/dense_matrix.h>

#include <math.h>
#include "opencl_adapter.hpp"

#define USE_DOUBLE_PRECISION
#define DEFAULT_TOLLERANCE 1e-12

// #ifndef m_kernel__
// #define m_kernel__
// #define m_global__
// #define m_local__
// #define m_constant__
// #endif

namespace libMesh {
    class Elem;
}

namespace utopia {

    class Intersector : public ::moonolith::OpenCLAdapter {
    public:
        // I do not know why the compiler wants this...
        template <typename T>
        inline static T sqrt(const T v) {
            return std::sqrt(v);
        }

#include "all_kernels.cl"
    };

    typedef utopia::Intersector::PMesh Polyhedron;

    class Transform;
    class QMortar;

    bool intersect_convex_polygons(const int n_vertices_1,
                                   const double *polygon_1,
                                   const int n_vertices_2,
                                   const double *polygon_2,
                                   int *n_vertices_result,
                                   double *result_buffer,
                                   double tol);

    void make_composite_quadrature_2D_non_affine(const libMesh::DenseMatrix<libMesh::Real> &polygon,
                                                 const double weight,
                                                 const int order,
                                                 QMortar &c_ir);
    void make_composite_quadrature_2D_from_tri_mesh(const std::vector<int> &tri,
                                                    const libMesh::DenseMatrix<libMesh::Real> &points,
                                                    const double weight,
                                                    const int order,
                                                    QMortar &c_ir);
    void make_composite_quadrature_2D(const libMesh::DenseMatrix<libMesh::Real> &polygon,
                                      const double weight,
                                      const int order,
                                      QMortar &c_ir);
    void make_composite_quadrature_3D(const Polyhedron &polyhedron,
                                      const double weight,
                                      const int order,
                                      QMortar &c_ir);
    void make_composite_quadrature_on_surf_2D(const libMesh::DenseMatrix<libMesh::Real> &line,
                                              const double weight,
                                              const int order,
                                              QMortar &c_ir);
    void make_composite_quadrature_on_surf_3D(const libMesh::DenseMatrix<libMesh::Real> &polygon,
                                              const double weight,
                                              const int order,
                                              QMortar &c_ir);

    void transform_to_reference(const Transform &trans, const int type, const QMortar &global_ir, QMortar &ref_ir);
    void transform_to_reference_surf(const Transform &trans, const int type, const QMortar &global_ir, QMortar &ref_ir);
    double compute_volume(const Polyhedron &poly);

    void make_polyhedron(const libMesh::Elem &e, Polyhedron &polyhedron);
    void make_polyline(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polyline);
    void make_polygon(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon);
    void make_polygon_3(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon);

    bool intersect_2D(const libMesh::DenseMatrix<libMesh::Real> &poly1,
                      const libMesh::DenseMatrix<libMesh::Real> &poly2,
                      libMesh::DenseMatrix<libMesh::Real> &intersection);
    bool intersect_3D(const libMesh::Elem &el1, const libMesh::Elem &el2, Polyhedron &intersection);
    bool intersect_3D(const Polyhedron &poly1, const Polyhedron &poly2, Polyhedron &intersection);

    bool project_2D(const libMesh::DenseMatrix<libMesh::Real> &poly1,
                    const libMesh::DenseMatrix<libMesh::Real> &poly2,
                    libMesh::DenseMatrix<libMesh::Real> &projection_1,
                    libMesh::DenseMatrix<libMesh::Real> &projection_2);

    bool project_3D(const libMesh::DenseMatrix<libMesh::Real> &polygon_1,
                    const libMesh::DenseMatrix<libMesh::Real> &polygon_2,
                    libMesh::DenseMatrix<libMesh::Real> &projection_1,
                    libMesh::DenseMatrix<libMesh::Real> &projection_2);

}  // namespace utopia

#undef mortar_assemble

#endif  // UTOPIA_INTERSECTOR_HPP

#ifndef UTOPIA_LIBMESH_DEPRECATED_HPP
#define UTOPIA_LIBMESH_DEPRECATED_HPP

#include "moonolith_function_space.hpp"
#include "utopia_libmesh_Transform.hpp"

#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_moonolith_libmesh_Convert.hpp"

#include "utopia_intersector.hpp"

namespace utopia {

    class LibMeshFunctionSpace;

    void compute_side_normal(const int dim, const libMesh::Elem &side, libMesh::Point &n);

    int order_for_l2_integral(const int dim,
                              const libMesh::Elem &master_el,
                              const int master_order,
                              const libMesh::Elem &slave_el,
                              const int slave_order);

    void print(const libMesh::QBase &ir, std::ostream &os = std::cout);
    double sum_of_weights(const libMesh::QBase &ir);
    double sum(const libMesh::DenseMatrix<libMesh::Real> &mat);

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

    void mortar_assemble(const libMesh::FEBase &trial_fe,
                         const libMesh::FEBase &test_fe,
                         libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble(const libMesh::FEVectorBase &trial_fe,
                         const libMesh::FEVectorBase &test_fe,
                         libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_normal_and_gap_assemble(const libMesh::FEBase &test_fe,
                                        const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                        const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble(const uint dim,
                                        const libMesh::FEBase &test_fe,
                                        const libMesh::Point &surf_normal,
                                        const libMesh::Point &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble(const libMesh::FEVectorBase &test_fe,
                                        const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                        const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble(const uint dim,
                                        const libMesh::FEVectorBase &test_fe,
                                        const libMesh::Point &surf_normal,
                                        const libMesh::Point &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const libMesh::FEBase &test_fe,
                                               const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                               const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const uint dim,
                                               const libMesh::FEBase &test_fe,
                                               const libMesh::Point &surf_normal,
                                               const libMesh::Point &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const libMesh::FEVectorBase &test_fe,
                                               const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                               const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const uint dim,
                                               const libMesh::FEVectorBase &test_fe,
                                               const libMesh::Point &surf_normal,
                                               const libMesh::Point &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap);

    // bool mortar_assemble(LibMeshFESpaceBase &src, LibMeshFESpaceBase &dest,
    // std::shared_ptr<libMesh::SparseMatrix<libMesh::Real> > &B);

    // bool transfer(LibMeshFESpaceBase &src, libMesh::DenseVector<libMesh::Real> &src_fun, LibMeshFESpaceBase &dest,
    // libMesh::DenseVector<libMesh::Real> &dest_fun);

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

    bool biorthgonal_weights(const int type, libMesh::Real &w_ii, libMesh::Real &w_ij);

    void mortar_assemble_biorth(const libMesh::FEBase &trial_fe,
                                const libMesh::FEBase &test_fe,
                                const int type,
                                libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble_biorth(const libMesh::FEVectorBase &trial_fe,
                                const libMesh::FEVectorBase &test_fe,
                                const int type,
                                libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble_biorth(const int dim,
                                const libMesh::FEBase &trial_fe,
                                const libMesh::FEBase &test_fe,
                                const int type,
                                const libMesh::DenseVector<libMesh::Real> &indicator,
                                libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble_biorth(const int dim,
                                const libMesh::FEVectorBase &trial_fe,
                                const libMesh::FEVectorBase &test_fe,
                                const int type,
                                const libMesh::DenseVector<libMesh::Real> &indicator,
                                libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble_weights(const libMesh::FEVectorBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights);
    void mortar_assemble_weights(const libMesh::FEBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights);

    void mortar_assemble_weighted_biorth(const libMesh::FEBase &trial_fe,
                                         const libMesh::FEBase &test_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                                         libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble_weighted_biorth(const libMesh::FEVectorBase &trial_fe,
                                         const libMesh::FEVectorBase &test_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                                         libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_normal_and_gap_assemble_weighted_biorth(const libMesh::FEVectorBase &test_fe,
                                                        const int dim,
                                                        const libMesh::Point &surf_normal,
                                                        const libMesh::Point &plane_normal,
                                                        const libMesh::Real &plane_offset,
                                                        const libMesh::DenseMatrix<libMesh::Real> &weights,
                                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                                        libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble_weighted_biorth(const libMesh::FEBase &test_fe,
                                                        const int dim,
                                                        const libMesh::Point &surf_normal,
                                                        const libMesh::Point &plane_normal,
                                                        const libMesh::Real &plane_offset,
                                                        const libMesh::DenseMatrix<libMesh::Real> &weights,
                                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                                        libMesh::DenseVector<libMesh::Real> &gap,
                                                        const bool visdebug = false);

    // void integrate_scalar_function(
    //     const libMesh::FEBase &test_fe,
    //     const std::vector<double> &fun,
    //     libMesh::DenseVector<libMesh::Real> &result
    // );

    template <class Array>
    void integrate_scalar_function(const libMesh::FEBase &test_fe,
                                   const Array &fun,
                                   libMesh::DenseVector<libMesh::Real> &result) {
        const auto &phi = test_fe.get_phi();
        const auto &dx = test_fe.get_JxW();
        const auto n_qp = fun.size();
        const auto n_shape_functions = phi.size();

        assert(n_qp == phi[0].size());
        assert(n_qp == dx.size());

        result.resize(n_shape_functions);
        result.zero();

        for (std::size_t i = 0; i < n_shape_functions; ++i) {
            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                result(i) += phi[i][qp] * fun[qp] * dx[qp];
            }
        }
    }

    void integrate_point_function(const int dim,
                                  const libMesh::FEBase &test_fe,
                                  const std::vector<libMesh::Point> &fun,
                                  libMesh::DenseMatrix<libMesh::Real> &result);

    template <int Dim>
    void l2_project_normal(const libMesh::FEBase &test_fe,
                           const ::moonolith::Vector<double, Dim> &normal,
                           libMesh::DenseVector<libMesh::Real> &result) {
        const auto &phi = test_fe.get_phi();
        const auto &dx = test_fe.get_JxW();
        const auto n_qp = dx.size();
        const auto n_shape_functions = phi.size();

        assert(n_qp == phi[0].size());

        result.resize(n_shape_functions * Dim);
        result.zero();

        for (std::size_t i = 0; i < n_shape_functions; ++i) {
            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                for (std::size_t d = 0; d < Dim; ++d) {
                    result(i * Dim + d) += phi[i][qp] * normal[d] * dx[qp];
                }
            }
        }
    }

    inline double biorth(const std::vector<std::vector<double>> &phi,
                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                         const int i,
                         const int qp) {
        double ret = 0.;
        auto n_funs = phi.size();
        for (std::size_t k = 0; k < n_funs; ++k) {
            ret += phi[k][qp] * weights(i, k);
        }

        return ret;
    }

    template <class Array>
    void integrate_scalar_function_weighted_biorth(const libMesh::FEBase &test_fe,
                                                   const libMesh::DenseMatrix<libMesh::Real> &weights,
                                                   const Array &fun,
                                                   libMesh::DenseVector<libMesh::Real> &result) {
        const auto &phi = test_fe.get_phi();
        const auto &dx = test_fe.get_JxW();
        const auto n_qp = fun.size();
        const auto n_shape_functions = phi.size();

        assert(n_qp == phi[0].size());
        assert(n_qp == dx.size());

        result.resize(n_shape_functions);
        result.zero();

        for (std::size_t i = 0; i < n_shape_functions; ++i) {
            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                auto f = biorth(phi, weights, i, qp);
                result(i) += f * fun[qp] * dx[qp];
            }
        }
    }

    template <int Dim>
    void l2_project_normal_weighted_biorth(const libMesh::FEBase &test_fe,
                                           const libMesh::DenseMatrix<libMesh::Real> &weights,
                                           const ::moonolith::Vector<double, Dim> &normal,

                                           libMesh::DenseVector<libMesh::Real> &result) {
        const auto &phi = test_fe.get_phi();
        const auto &dx = test_fe.get_JxW();
        const auto n_qp = dx.size();
        const auto n_shape_functions = phi.size();

        assert(n_qp == phi[0].size());

        result.resize(n_shape_functions * Dim);
        result.zero();

        for (std::size_t i = 0; i < n_shape_functions; ++i) {
            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                auto f = biorth(phi, weights, i, qp);
                for (std::size_t d = 0; d < Dim; ++d) {
                    result(i * Dim + d) += f * normal[d] * dx[qp];
                }
            }
        }
    }

    void mortar_assemble_weighted_biorth(const libMesh::FEBase &trial_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &trafo,
                                         const libMesh::FEBase &test_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                                         libMesh::DenseMatrix<libMesh::Real> &elmat);

    // template <int Dim>
    // class ConvertFunctionSpace<LibMeshFunctionSpace, ::moonolith::FunctionSpace<::moonolith::Mesh<double, Dim>>> {
    // public:
    //     inline static void apply(const LibMeshFunctionSpace &in,
    //                              ::moonolith::FunctionSpace<::moonolith::Mesh<double, Dim>> &out) {
    //         int var_num = in.subspace_id();
    //         convert_libmesh_to_moonolith(in.mesh(), in.dof_map(), var_num, out);
    //     }

    //     inline static void apply(const libMesh::MeshBase &in,
    //                              const libMesh::DofMap &dof_map,
    //                              unsigned int var_num,
    //                              ::moonolith::FunctionSpace<::moonolith::Mesh<double, Dim>> &out) {
    //         convert_libmesh_to_moonolith(in, dof_map, var_num, out);
    //     }
    // };

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_DEPRECATED_HPP

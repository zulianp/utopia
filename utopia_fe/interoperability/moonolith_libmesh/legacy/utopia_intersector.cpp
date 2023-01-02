
#include <iostream>
#include "clipper.hpp"

#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/quadrature_gauss.h>

#include "Box.hpp"
#include "utopia_libmesh_QMortar.hpp"

#include "utopia_Polygon.hpp"

#include "utopia_intersector.hpp"

// #define USE_CLIPPER 1

namespace utopia {

    template <class V, class T>
    static void add(const V &p, const T &alpha, const V &v, V &result) {
        assert(result.size() == p.size());

        for (int i = 0; i < result.size(); ++i) {
            result(i) = p(i) + alpha * v(i);
        }
    }

    bool intersect_convex_polygons(const int n_vertices_1,
                                   const double *polygon_1,
                                   const int n_vertices_2,
                                   const double *polygon_2,
                                   int *n_vertices_result,
                                   double *result_buffer,
                                   double tol) {
        // static bool lib_msg = false;

#ifdef USE_CLIPPER
        // if(!lib_msg) {
        // 	std::cout << "using clipper" << std::endl;
        // 	lib_msg = true;
        // }

        using namespace ClipperLib;

        double min_x = polygon_1[0];
        double min_y = polygon_1[1];
        double max_x = polygon_1[0];
        double max_y = polygon_1[1];

        for (int i = 1; i < n_vertices_1; ++i) {
            const int i2 = i * 2;
            min_x = std::min(min_x, polygon_1[i2]);
            max_x = std::max(max_x, polygon_1[i2]);

            min_y = std::min(min_y, polygon_1[i2 + 1]);
            max_y = std::max(max_y, polygon_1[i2 + 1]);
        }

        for (int i = 0; i < n_vertices_2; ++i) {
            const int i2 = i * 2;
            min_x = std::min(min_x, polygon_2[i2]);
            max_x = std::max(max_x, polygon_2[i2]);

            min_y = std::min(min_y, polygon_2[i2 + 1]);
            max_y = std::max(max_y, polygon_2[i2 + 1]);
        }

        const double cut_off = 1e10 / std::max(max_x - min_x, max_y - min_y);  // 1e18 is the max represented value 

        Paths subj(1), clip(1), solution;

        subj[0].reserve(n_vertices_1);
        clip[0].reserve(n_vertices_2);

        for (int i = 0; i < n_vertices_1; ++i) {
            const int i2 = i * 2;
            subj[0].push_back(IntPoint((polygon_1[i2] - min_x) * cut_off, (polygon_1[i2 + 1] - min_y) * cut_off));
        }

        for (int i = 0; i < n_vertices_2; ++i) {
            const int i2 = i * 2;
            clip[0].push_back(IntPoint((polygon_2[i2] - min_x) * cut_off, (polygon_2[i2 + 1] - min_y) * cut_off));
        }

        // perform intersection ...
        Clipper c;
        c.AddPaths(subj, ptSubject, true);
        c.AddPaths(clip, ptClip, true);
        c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);

        if (solution.empty() || solution[0].size() < 3) return false;

        assert(solution.size() == 1);

        *n_vertices_result = solution[0].size();

        for (std::size_t i = 0; i < solution[0].size(); ++i) {
            const int i2 = i * 2;
            result_buffer[i2] = solution[0][i].X / cut_off + min_x;
            result_buffer[i2 + 1] = solution[0][i].Y / cut_off + min_y;
        }

        return true;

#else   // USE_CLIPPER

        // if(!lib_msg) {
        // 	std::cout << "using homemade" << std::endl;
        // 	lib_msg = true;
        // }

        return Intersector::intersect_convex_polygons(
            n_vertices_1, polygon_1, n_vertices_2, polygon_2, n_vertices_result, result_buffer, tol);
#endif  // USE_CLIPPER
    }

    void make_composite_quadrature_2D_non_affine(const libMesh::DenseMatrix<libMesh::Real> &polygon,
                                                 const double weight,
                                                 const int order,
                                                 QMortar &c_ir) {
        std::vector<int> tri;
        // triangulate_polygon(polygon.m(), &polygon.get_values()[0], tri);
        // make_composite_quadrature_2D_from_tri_mesh(tri, polygon, weight, order, c_ir);
        assert(false);
    }

    void make_composite_quadrature_2D_from_tri_mesh(const std::vector<int> &tri,
                                                    const libMesh::DenseMatrix<libMesh::Real> &points,
                                                    const double weight,
                                                    const int order,
                                                    QMortar &c_ir) {
        using namespace libMesh;
        // FIXME only for sligthlty deformed elements
        libMesh::QGauss ir(2, libMesh::Order(order));
        ir.init(libMesh::TRI6);

        double normalization_factor = weight * 2.;

        double triangle[3 * 2] = {0., 0., 0., 0., 0., 0.};

        const int n_triangles = tri.size() / 3;
        const int n_quad_points = n_triangles * ir.n_points();
        c_ir.resize(n_quad_points);

        libMesh::DenseVector<libMesh::Real> o, u, v, p;
        double relative_weight = 0;
        int quad_index = 0;
        for (int i = 0; i < n_triangles; ++i) {
            const int i3 = i * 3;
            int v1 = tri[i3];
            int v2 = tri[i3 + 1];
            int v3 = tri[i3 + 2];

            get_row(v1, points, o);
            get_row(v2, points, u);
            get_row(v3, points, v);

            triangle[0] = o(0);
            triangle[1] = o(1);

            triangle[2] = u(0);
            triangle[3] = u(1);

            triangle[4] = v(0);
            triangle[5] = v(1);

            // plot_polygon(2, 3, triangle);

            u -= o;
            v -= o;

            const double scale = fabs(Intersector::polygon_area_2(3, triangle)) * normalization_factor;

            relative_weight += scale;

            for (int k = 0; k < ir.n_points(); ++k, ++quad_index) {
                auto &qp = ir.get_points()[k];
                auto &c_qp = c_ir.get_points()[quad_index];

                p = o;
                add(p, qp(0), u, p);
                add(p, qp(1), v, p);

                c_qp(0) = p(0);
                c_qp(1) = p(1);
                c_qp(2) = 0.0;

                c_ir.get_weights()[quad_index] = ir.w(k) * scale;
            }
        }

        assert(relative_weight <= 1.0001);
        assert(quad_index == n_quad_points);
        // plot_quad_points(2, c_ir.get_points(), "qp");
    }

    void make_composite_quadrature_2D(const libMesh::DenseMatrix<libMesh::Real> &polygon,
                                      const double weight,
                                      const int order,
                                      QMortar &c_ir) {
        libMesh::QGauss ir(2, libMesh::Order(order));
        ir.init(libMesh::TRI3);

        const double normalization_factor = weight * 2.;

        const int n_triangles = polygon.m() - 2;
        const int n_quad_points = n_triangles * ir.n_points();

        assert(fabs(sum_of_weights(ir) - 0.5) < 1e-8);

        double triangle[3 * 2] = {polygon.get_values()[0], polygon.get_values()[1], 0., 0., 0., 0.};

        libMesh::DenseVector<libMesh::Real> o, u, v, p;
        get_row(0, polygon, o);

        c_ir.resize(n_quad_points);

        double sum_weights = 0.;

        int quad_index = 0;
        for (int i = 2; i < polygon.m(); ++i) {
            get_row(i - 1, polygon, u);
            get_row(i, polygon, v);

            triangle[2] = u(0);
            triangle[3] = u(1);

            triangle[4] = v(0);
            triangle[5] = v(1);

            u -= o;
            v -= o;

            const double tri_area = Intersector::polygon_area_2(3, triangle);
            const double scale = std::abs(tri_area) * normalization_factor;

            assert(tri_area > 0.);
            assert(tri_area * normalization_factor <= 2.0001);

            for (int k = 0; k < ir.n_points(); ++k, ++quad_index) {
                auto &qp = ir.get_points()[k];
                auto &c_qp = c_ir.get_points()[quad_index];

                p = o;
                add(p, qp(0), u, p);
                add(p, qp(1), v, p);

                c_qp(0) = p(0);
                c_qp(1) = p(1);
                c_qp(2) = 0.0;

                c_ir.get_weights()[quad_index] = ir.w(k) * scale;
                sum_weights += c_ir.get_weights()[quad_index];
            }
        }

        assert(sum_weights > 0.);
        assert(sum_weights <= 1.0001);
        assert(quad_index == n_quad_points);
    }

    void make_composite_quadrature_3D(const Polyhedron &polyhedron,
                                      const double weight,
                                      const int order,
                                      QMortar &c_ir) {
        using std::min;

        libMesh::QGauss ir(3, libMesh::Order(order));
        ir.init(libMesh::TET4);

        const double normalization_factor = weight * 6.;

        const int n_sub_elements = Intersector::n_volume_elements(polyhedron);
        const int total_n_quad_points = n_sub_elements * ir.n_points();

        std::vector<double> ref_quad_points(ir.n_points() * 3);
        std::vector<double> ref_quad_weights(ir.n_points());

        for (int i = 0; i < ir.n_points(); ++i) {
            const int offset = i * 3;

            ref_quad_points[offset] = ir.qp(i)(0);
            ref_quad_points[offset + 1] = ir.qp(i)(1);
            ref_quad_points[offset + 2] = ir.qp(i)(2);

            ref_quad_weights[i] = ir.w(i);
        }

        c_ir.resize(total_n_quad_points);

        // global quadrature points
        double quad_points[MAX_QUAD_POINTS * 3];
        double quad_weights[MAX_QUAD_POINTS];

        double barycenter_p[3];

        Intersector::row_average(polyhedron.n_nodes, polyhedron.n_dims, polyhedron.points, barycenter_p);

        const uint max_n_sub_els = Intersector::max_n_elements_from_facets(polyhedron);
        const uint max_n_sub_inc = std::max(1u, MAX_QUAD_POINTS / (max_n_sub_els * ir.n_points()));

        int utopia_fe_quad_index = 0;

        if (n_sub_elements == 1) {
            Intersector::tetrahedron_transform(
                polyhedron.points, total_n_quad_points, &ref_quad_points[0], quad_points);
            const double w = fabs(Intersector::m_tetrahedron_volume(polyhedron.points) * normalization_factor);

            for (int i = 0; i < total_n_quad_points; ++i, ++utopia_fe_quad_index) {
                const int offset = i * 3;
                c_ir.get_points()[utopia_fe_quad_index](0) = quad_points[offset];
                c_ir.get_points()[utopia_fe_quad_index](1) = quad_points[offset + 1];
                c_ir.get_points()[utopia_fe_quad_index](2) = quad_points[offset + 2];

                c_ir.get_weights()[utopia_fe_quad_index] = ref_quad_weights[i] * w;
            }

        } else {
            for (int begin_k = 0; begin_k < polyhedron.n_elements;) {
                const int end_k = min(begin_k + max_n_sub_inc, static_cast<uint>(polyhedron.n_elements));
                assert(end_k > begin_k && "end_k > begin_k");

                const int n_quad_points =
                    Intersector::make_quadrature_points_from_polyhedron_in_range_around_point(polyhedron,
                                                                                              begin_k,
                                                                                              end_k,
                                                                                              normalization_factor,
                                                                                              ir.n_points(),
                                                                                              &ref_quad_points[0],
                                                                                              &ref_quad_weights[0],
                                                                                              barycenter_p,
                                                                                              quad_points,
                                                                                              quad_weights);

                for (int i = 0; i < n_quad_points; ++i, ++utopia_fe_quad_index) {
                    const int offset = i * 3;
                    c_ir.get_points()[utopia_fe_quad_index](0) = quad_points[offset];
                    c_ir.get_points()[utopia_fe_quad_index](1) = quad_points[offset + 1];
                    c_ir.get_points()[utopia_fe_quad_index](2) = quad_points[offset + 2];

                    c_ir.get_weights()[utopia_fe_quad_index] = quad_weights[i];
                }

                begin_k = end_k;
            }
        }
    }

    void make_composite_quadrature_on_surf_2D(const libMesh::DenseMatrix<libMesh::Real> &line,
                                              const double weight,
                                              const int order,
                                              QMortar &c_ir) {
        using namespace libMesh;
        QGauss ir(1, libMesh::Order(order));
        ir.init(libMesh::EDGE2);

        c_ir.resize(ir.n_points());

        const double normalization_factor = weight * 0.5;

        Point o, v;

        for (int i = 0; i < line.n(); ++i) {
            o(i) = line(0, i);
            v(i) = line(1, i) - line(0, i);
        }

        const Real length = v.norm();

        for (int k = 0; k < ir.n_points(); ++k) {
            const Real w = (1 + ir.get_points()[k](0)) * 0.5;
            assert(w >= -1e-16);
            assert(w <= 1 + 1e-16);

            c_ir.get_points()[k] = v;
            c_ir.get_points()[k] *= w;
            c_ir.get_points()[k] += o;

            c_ir.get_weights()[k] = ir.w(k) * length * normalization_factor;
        }
    }

    void make_composite_quadrature_on_surf_3D(const libMesh::DenseMatrix<libMesh::Real> &polygon,
                                              const double weight,
                                              const int order,
                                              QMortar &c_ir) {
        libMesh::QGauss ir(2, libMesh::Order(order));

        if (order <= 2) {
            ir.init(libMesh::TRI3);
        } else {
            ir.init(libMesh::TRI6);
        }

        // remove scaling of triangle quad rule..
        const double normalization_factor = weight * 2.;

        const int n_triangles = polygon.m() - 2;
        const int n_quad_points = n_triangles * ir.n_points();

        assert(fabs(sum_of_weights(ir) - 0.5) < 1e-8);

        double triangle[3 * 3] = {
            polygon.get_values()[0], polygon.get_values()[1], polygon.get_values()[2], 0., 0., 0., 0., 0., 0.};

        libMesh::DenseVector<libMesh::Real> o, u, v, p;
        get_row(0, polygon, o);

        c_ir.resize(n_quad_points);

        double relative_weight = 0;

        int quad_index = 0;
        for (int i = 2; i < polygon.m(); ++i) {
            get_row(i - 1, polygon, u);
            get_row(i, polygon, v);

            triangle[3] = u(0);
            triangle[4] = u(1);
            triangle[5] = u(2);

            triangle[6] = v(0);
            triangle[7] = v(1);
            triangle[8] = v(2);

            u -= o;
            v -= o;

            relative_weight = Intersector::polygon_area_3(3, triangle) * normalization_factor;

            assert(relative_weight >= 0.);

            for (int k = 0; k < ir.n_points(); ++k, ++quad_index) {
                auto &qp = ir.get_points()[k];
                auto &c_qp = c_ir.get_points()[quad_index];

                p = o;
                add(p, qp(0), u, p);
                add(p, qp(1), v, p);

                c_qp(0) = p(0);
                c_qp(1) = p(1);
                c_qp(2) = p(2);

                c_ir.get_weights()[quad_index] = ir.w(k) * relative_weight;
            }
        }

        assert(quad_index == n_quad_points);
    }

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const libMesh::FEBase &test_fe,
                                               const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                               const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap) {
        using namespace libMesh;

        const uint dim = plane_normal.size();

        libMesh::Real w_ii, w_ij;
        biorthgonal_weights(type, w_ii, w_ij);

        if (normals.m() != test_fe.get_phi().size() || dim != normals.n()) {
            normals.resize(test_fe.get_phi().size(), dim);
            normals.zero();
            gap.resize(test_fe.get_phi().size());
            gap.zero();
        }

        const auto &test = test_fe.get_phi();
        // const auto &grad   = test_fe.get_dphi();
        const auto &point = test_fe.get_xyz();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_qp = test[0].size();

        DenseVector<Real> p(dim);
        // DenseVector<Real> v(dim);

        for (uint i = 0; i < n_test; ++i) {
            for (uint qp = 0; qp < n_qp; ++qp) {
                p(0) = point[qp](0);
                p(1) = point[qp](1);

                if (dim > 2) {
                    p(2) = point[qp](2);
                }

                Real isect = 0;

                Intersector::intersect_ray_with_plane(dim,
                                                      1,
                                                      &p.get_values()[0],
                                                      &surf_normal.get_values()[0],
                                                      &plane_normal.get_values()[0],
                                                      plane_offset,
                                                      &isect);

                // v = surf_normal;
                // v *= isect;
                // quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]);

                // printf("g: %g (%g, %g)\n", isect, p(1), plane_offset);
                // assert(isect > 0);

                auto biorth_test = ((0 == i) ? w_ii : w_ij) * test[0][qp];

                for (uint k = 1; k < n_test; ++k) {
                    biorth_test += ((k == i) ? w_ii : w_ij) * test[k][qp];
                }

                for (int k = 0; k < dim; ++k) {
                    if (i % dim == k) {
                        make_tp(k, biorth_test);
                        break;
                    }
                }

                gap(i) += indicator(i) * biorth_test * isect * JxW[qp];

                for (uint d = 0; d < dim; ++d) {
                    normals(i, d) += indicator(i) * biorth_test * surf_normal(d) * JxW[qp];
                }
            }
        }
    }

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const libMesh::FEVectorBase &test_fe,
                                               const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                               const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap) {
        libMesh::Real w_ii, w_ij;
        biorthgonal_weights(type, w_ii, w_ij);

        // std::cout << w_ii << " " << w_ij << std::endl;

        using namespace libMesh;

        const uint dim = plane_normal.size();

        if (normals.m() != test_fe.get_phi().size() / dim || dim != normals.n()) {
            normals.resize(test_fe.get_phi().size() / dim, dim);
            normals.zero();
            gap.resize(test_fe.get_phi().size());
            gap.zero();
        }

        const auto &test = test_fe.get_phi();
        const auto &point = test_fe.get_xyz();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_qp = test[0].size();

        DenseVector<Real> p(dim);
        DenseVector<Real> v(dim);

        for (uint qp = 0; qp < n_qp; ++qp) {
            p(0) = point[qp](0);
            p(1) = point[qp](1);

            if (dim > 2) {
                p(2) = point[qp](2);
            }

            Real isect = 0;

            Intersector::intersect_ray_with_plane(dim,
                                                  1,
                                                  &p.get_values()[0],
                                                  &surf_normal.get_values()[0],
                                                  &plane_normal.get_values()[0],
                                                  plane_offset,
                                                  &isect);

            v = surf_normal;
            v *= isect;
            // quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]);

            for (uint i = 0; i < n_test; ++i) {
                auto biorth_test = ((0 == i) ? w_ii : w_ij) * test[0][qp];

                for (uint k = 1; k < n_test; ++k) {
                    biorth_test += ((k == i) ? w_ii : w_ij) * test[k][qp];
                }

                for (int k = 0; k < dim; ++k) {
                    if (i % dim == k) {
                        make_tp(k, biorth_test);
                        break;
                    }
                }

                gap(i) += indicator(i) * biorth_test(0) * isect * JxW[qp];

                for (uint d = 0; d < dim; ++d) {
                    normals.get_values()[i] += indicator(i) * biorth_test(d) * surf_normal(d) * JxW[qp];
                }
            }
        }
    }

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const uint dim,
                                               const libMesh::FEVectorBase &test_fe,
                                               const libMesh::Point &surf_normal,
                                               const libMesh::Point &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap) {
        using namespace libMesh;
        DenseVector<Real> surf_normal_v(dim), plane_normal_v(dim);

        for (uint i = 0; i < dim; ++i) {
            surf_normal_v(i) = surf_normal(i);
            plane_normal_v(i) = plane_normal(i);
        }

        mortar_normal_and_gap_assemble_biorth(
            type, test_fe, surf_normal_v, plane_normal_v, plane_offset, indicator, normals, gap);
    }

    bool intersect_2D(const libMesh::DenseMatrix<libMesh::Real> &poly1,
                      const libMesh::DenseMatrix<libMesh::Real> &poly2,
                      libMesh::DenseMatrix<libMesh::Real> &intersection) {
        double result_buffer[MAX_N_ISECT_POINTS * 2];
        int n_vertices_result;

        assert(Intersector::polygon_area_2(poly1.m(), &poly1.get_values()[0]) > 0);
        assert(Intersector::polygon_area_2(poly2.m(), &poly2.get_values()[0]) > 0);

        if (!
            // isector.
            intersect_convex_polygons(poly1.m(),
                                      &poly1.get_values()[0],
                                      poly2.m(),
                                      &poly2.get_values()[0],
                                      &n_vertices_result,
                                      result_buffer,
                                      DEFAULT_TOLLERANCE)) {
            return false;
        }

        assert(Intersector::polygon_area_2(n_vertices_result, result_buffer) > 0);

        intersection.resize(n_vertices_result, 2);
        std::copy(result_buffer, result_buffer + n_vertices_result * 2, &intersection.get_values()[0]);

        // plot_polygon(2, n_vertices_result, &intersection.get_values()[0], "isect/poly");
        return true;
    }

    double compute_volume(const Polyhedron &poly) {
        if (poly.type >= P_MESH_TYPE_SURF) {
            return Intersector::polygon_area_3(poly.n_nodes, poly.points);
        } else {
            return Intersector::p_mesh_volume_3(poly);
        }
    }

    bool intersect_3D(const Polyhedron &poly1, const Polyhedron &poly2, Polyhedron &intersection) {
        if (poly1.type >= P_MESH_TYPE_SURF) {
            bool ok = Intersector::intersect_convex_polyhedron_with_polygon(
                poly2, poly1.n_nodes, poly1.points, &intersection);
            intersection.type = P_MESH_TYPE_SURF;
            return ok;
        } else if (poly2.type >= P_MESH_TYPE_SURF) {
            bool ok = Intersector::intersect_convex_polyhedron_with_polygon(
                poly1, poly2.n_nodes, poly2.points, &intersection);
            intersection.type = P_MESH_TYPE_SURF;
            return ok;
        } else {
            bool ok = Intersector::intersect_convex_polyhedra(poly1, poly2, &intersection);
            intersection.type = P_MESH_TYPE_UNSTRUCTURED;
            return ok;
        }
    }

    bool intersect_3D(const libMesh::Elem &el1, const libMesh::Elem &el2, Polyhedron &intersection) {
        Polyhedron p1, p2;
        make_polyhedron(el1, p1);
        make_polyhedron(el2, p2);
        return intersect_3D(p1, p2, intersection);
    }

    static void make_polyhedron_from_generic_tet(const libMesh::Elem &e, Polyhedron &polyhedron) {
        polyhedron.n_elements = 4;
        polyhedron.n_nodes = 4;
        polyhedron.n_dims = 3;

        for (int i = 0; i < 4; ++i) {
            const int offset = i * 3;

            for (int j = 0; j < 3; ++j) {
                polyhedron.points[offset + j] = e.point(i)(j);
            }
        }

        polyhedron.el_ptr[0] = 0;
        polyhedron.el_ptr[1] = 3;
        polyhedron.el_ptr[2] = 6;
        polyhedron.el_ptr[3] = 9;
        polyhedron.el_ptr[4] = 12;

        // face 0
        polyhedron.el_index[0] = 0;
        polyhedron.el_index[1] = 1;
        polyhedron.el_index[2] = 3;

        // face 1
        polyhedron.el_index[3] = 1;
        polyhedron.el_index[4] = 2;
        polyhedron.el_index[5] = 3;

        // face 2
        polyhedron.el_index[6] = 0;
        polyhedron.el_index[7] = 3;
        polyhedron.el_index[8] = 2;

        // face 3
        polyhedron.el_index[9] = 1;
        polyhedron.el_index[10] = 2;
        polyhedron.el_index[11] = 0;

        polyhedron.type = P_MESH_TYPE_TET;
    }

    void make_polyhedron_from_tet4(const libMesh::Elem &e, Polyhedron &polyhedron) {
        assert(e.dim() == 3);
        assert(e.n_nodes() == 4);

        make_polyhedron_from_generic_tet(e, polyhedron);
    }

    void make_polyhedron_from_tet10(const libMesh::Elem &e, Polyhedron &polyhedron) {
        assert(e.dim() == 3);
        assert(e.n_nodes() == 10);

        make_polyhedron_from_generic_tet(e, polyhedron);
    }

    void make_polyhedron_from_pyramid(const libMesh::Elem &e, Polyhedron &polyhedron) {
        assert(e.dim() == 3);
        assert(is_pyramid(e.type()));

        polyhedron.n_elements = 5;
        polyhedron.n_nodes = 5;
        polyhedron.n_dims = 3;

        for (int i = 0; i < 5; ++i) {
            const int offset = i * 3;

            for (int j = 0; j < 3; ++j) {
                polyhedron.points[offset + j] = e.point(i)(j);
            }
        }

        polyhedron.el_ptr[0] = 0;
        polyhedron.el_ptr[1] = 3;
        polyhedron.el_ptr[2] = 6;
        polyhedron.el_ptr[3] = 9;
        polyhedron.el_ptr[4] = 12;
        polyhedron.el_ptr[5] = 16;

        // face 0
        polyhedron.el_index[0] = 0;
        polyhedron.el_index[1] = 1;
        polyhedron.el_index[2] = 4;

        // face 1
        polyhedron.el_index[3] = 1;
        polyhedron.el_index[4] = 2;
        polyhedron.el_index[5] = 4;

        // face 2
        polyhedron.el_index[6] = 3;
        polyhedron.el_index[7] = 0;
        polyhedron.el_index[8] = 4;

        // face 3
        polyhedron.el_index[9] = 2;
        polyhedron.el_index[10] = 3;
        polyhedron.el_index[11] = 4;

        // face 4
        polyhedron.el_index[12] = 0;
        polyhedron.el_index[13] = 3;
        polyhedron.el_index[14] = 2;
        polyhedron.el_index[15] = 1;

        polyhedron.type = P_MESH_TYPE_UNSTRUCTURED;
    }

    void make_polyhedron_from_prism(const libMesh::Elem &e, Polyhedron &polyhedron) {
        assert(e.dim() == 3);
        assert(is_prism(e.type()));

        polyhedron.n_elements = 5;
        polyhedron.n_nodes = 6;
        polyhedron.n_dims = 3;

        for (int i = 0; i < 6; ++i) {
            const int offset = i * 3;

            for (int j = 0; j < 3; ++j) {
                polyhedron.points[offset + j] = e.point(i)(j);
            }
        }

        polyhedron.el_ptr[0] = 0;
        polyhedron.el_ptr[1] = 3;
        polyhedron.el_ptr[2] = 6;
        polyhedron.el_ptr[3] = 10;
        polyhedron.el_ptr[4] = 14;
        polyhedron.el_ptr[5] = 18;

        // face 0
        polyhedron.el_index[0] = 0;
        polyhedron.el_index[1] = 1;
        polyhedron.el_index[2] = 2;

        // face 1
        polyhedron.el_index[3] = 4;
        polyhedron.el_index[4] = 3;
        polyhedron.el_index[5] = 5;

        // face 2
        polyhedron.el_index[6] = 0;
        polyhedron.el_index[7] = 2;
        polyhedron.el_index[8] = 5;
        polyhedron.el_index[9] = 3;

        // face 3
        polyhedron.el_index[10] = 1;
        polyhedron.el_index[11] = 4;
        polyhedron.el_index[12] = 5;
        polyhedron.el_index[13] = 2;

        // face 4
        polyhedron.el_index[14] = 0;
        polyhedron.el_index[15] = 1;
        polyhedron.el_index[16] = 4;
        polyhedron.el_index[17] = 3;

        polyhedron.type = P_MESH_TYPE_UNSTRUCTURED;
    }

    static void make_polyhedron_from_shell_element(const libMesh::Elem &e, Polyhedron &polyhedron) {
        polyhedron.n_dims = 3;

        polyhedron.el_ptr[0] = 0;
        polyhedron.el_ptr[1] = 2;
        polyhedron.el_ptr[2] = 4;

        polyhedron.el_index[0] = 0;
        polyhedron.el_index[1] = 1;

        polyhedron.el_index[2] = 1;
        polyhedron.el_index[3] = 2;

        polyhedron.el_index[4] = 2;
        polyhedron.el_index[5] = 3;

        if (is_tri(e.type())) {
            polyhedron.n_elements = 3;
            polyhedron.n_nodes = 3;
            polyhedron.type = P_MESH_TYPE_TRIANGLE;

        } else if (is_quad(e.type())) {
            polyhedron.n_elements = 4;
            polyhedron.n_nodes = 4;
            polyhedron.el_ptr[3] = 6;

            polyhedron.el_index[6] = 3;
            polyhedron.el_index[7] = 4;

            polyhedron.type = P_MESH_TYPE_QUAD;

        } else {
            assert(false);
        }

        for (int i = 0; i < polyhedron.n_nodes; ++i) {
            const int offset = i * 3;

            for (int j = 0; j < 3; ++j) {
                polyhedron.points[offset + j] = e.point(i)(j);
            }
        }
    }

    static void make_polyhedron_from_generic_hex(const libMesh::Elem &e, Polyhedron &polyhedron) {
        polyhedron.n_elements = 6;
        polyhedron.n_nodes = 8;
        polyhedron.n_dims = 3;

        for (int i = 0; i < 8; ++i) {
            const int offset = i * 3;

            for (int j = 0; j < 3; ++j) {
                polyhedron.points[offset + j] = e.point(i)(j);
            }
        }

        polyhedron.el_ptr[0] = 0;
        polyhedron.el_ptr[1] = 4;
        polyhedron.el_ptr[2] = 8;
        polyhedron.el_ptr[3] = 12;
        polyhedron.el_ptr[4] = 16;
        polyhedron.el_ptr[5] = 20;
        polyhedron.el_ptr[6] = 24;

        // face 0
        polyhedron.el_index[0] = 0;
        polyhedron.el_index[1] = 1;
        polyhedron.el_index[2] = 5;
        polyhedron.el_index[3] = 4;

        // face 1
        polyhedron.el_index[4] = 1;
        polyhedron.el_index[5] = 2;
        polyhedron.el_index[6] = 6;
        polyhedron.el_index[7] = 5;

        // face 2
        polyhedron.el_index[8] = 3;
        polyhedron.el_index[9] = 7;
        polyhedron.el_index[10] = 6;
        polyhedron.el_index[11] = 2;

        // face 3
        polyhedron.el_index[12] = 0;
        polyhedron.el_index[13] = 4;
        polyhedron.el_index[14] = 7;
        polyhedron.el_index[15] = 3;

        // face 4
        polyhedron.el_index[16] = 2;
        polyhedron.el_index[17] = 1;
        polyhedron.el_index[18] = 0;
        polyhedron.el_index[19] = 3;

        // face 5
        polyhedron.el_index[20] = 6;
        polyhedron.el_index[21] = 7;
        polyhedron.el_index[22] = 4;
        polyhedron.el_index[23] = 5;

        polyhedron.type = P_MESH_TYPE_HEX;
    }

    // FIXME between -1, 1 ?
    void make_polyhedron_from_hex8(const libMesh::Elem &e, Polyhedron &polyhedron) {
        assert(e.dim() == 3);
        assert(e.n_nodes() == 8);

        make_polyhedron_from_generic_hex(e, polyhedron);
    }

    // FIXME between -1, 1 ?
    void make_polyhedron_from_hex27(const libMesh::Elem &e, Polyhedron &polyhedron) {
        assert(e.dim() == 3);
        assert(e.n_nodes() == 27);

        make_polyhedron_from_generic_hex(e, polyhedron);
    }

    void make_polyhedron(const libMesh::Elem &e, Polyhedron &polyhedron) {
        if (is_tri(e.type()) || is_quad(e.type())) {
            make_polyhedron_from_shell_element(e, polyhedron);
            return;
        }

        if (is_pyramid(e.type())) {
            make_polyhedron_from_pyramid(e, polyhedron);
            return;
        }

        if (is_prism(e.type())) {
            make_polyhedron_from_prism(e, polyhedron);
            return;
        }

        // FIXME use libMesh enum types
        switch (e.n_nodes()) {
            case 4: {
                make_polyhedron_from_tet4(e, polyhedron);
                break;
            }

            case 8: {
                make_polyhedron_from_hex8(e, polyhedron);
                break;
            }

            case 10: {
                make_polyhedron_from_tet10(e, polyhedron);
                break;
            }

            case 20: {
                make_polyhedron_from_generic_hex(e, polyhedron);
                break;
            }

            case 27: {
                make_polyhedron_from_hex27(e, polyhedron);
                break;
            }

            default: {
                assert(false);
                break;
            }
        }
    }

    bool project_2D(const libMesh::DenseMatrix<libMesh::Real> &poly1,
                    const libMesh::DenseMatrix<libMesh::Real> &poly2,
                    libMesh::DenseMatrix<libMesh::Real> &projection_1,
                    libMesh::DenseMatrix<libMesh::Real> &projection_2) {
        using std::max;
        using std::min;

        typedef Intersector::Scalar Scalar;

        Scalar A[2 * 2], b[2];
        Scalar Ainv[2 * 2], binv[2];

        //      Scalar normal_master[2];
        //      Scalar normal_slave [2];

        Intersector::SurfaceMortarWorkspace w;

        const uint n_points_master = poly1.m();

        //////////////////////////////////////////////////////////////////////////////////////////
        // computing geometric surface projection

        // moving from global space to reference space

        Intersector::line_make_affine_transform_2(&poly2.get_values()[0], A, b);
        Intersector::make_inverse_affine_transform_2(A, b, Ainv, binv);

        Intersector::apply_affine_transform_2(Ainv, binv, n_points_master, &poly1.get_values()[0], w.ref_points_master);

        Scalar x_min, y_min;
        Scalar x_max, y_max;

        if (w.ref_points_master[0] < w.ref_points_master[2]) {
            x_min = w.ref_points_master[0];
            y_min = w.ref_points_master[1];

            x_max = w.ref_points_master[2];
            y_max = w.ref_points_master[3];
        } else {
            x_min = w.ref_points_master[2];
            y_min = w.ref_points_master[3];

            x_max = w.ref_points_master[0];
            y_max = w.ref_points_master[1];
        }

        // check if there is any intersection
        if (x_max <= 0) {
            return false;
        }

        if (x_min >= 1) {
            return false;
        }

        const Scalar x_min_isect = max(x_min, (Scalar)0);
        const Scalar x_max_isect = min(x_max, (Scalar)1);

        const Scalar dx = (x_max - x_min);
        const Scalar dy = (y_max - y_min);

        const Scalar y_min_isect = (x_min_isect - x_min) / dx * dy + y_min;
        const Scalar y_max_isect = (x_max_isect - x_min) / dx * dy + y_min;

        // store intersection lines
        w.isect_master[0] = x_min_isect;
        w.isect_master[1] = y_min_isect;
        w.isect_master[2] = x_max_isect;
        w.isect_master[3] = y_max_isect;

        w.isect_slave[0] = x_min_isect;
        w.isect_slave[1] = 0;
        w.isect_slave[2] = x_max_isect;
        w.isect_slave[3] = 0;

        //      const Scalar inv_area_slave = 2.0/( Intersector::det_2(A) );

        //////////////////////////////////////////////////////////////////////////////////////////
        // create master fe object from intersection

        projection_1.resize(2, 2);
        projection_2.resize(2, 2);

        // move back to global coordinates

        Intersector::apply_affine_transform_2(A, b, 2, w.isect_master, &projection_1.get_values()[0]);

        // move back to global coordinates
        Intersector::apply_affine_transform_2(A, b, 2, w.isect_slave, &projection_2.get_values()[0]);
        return true;
    }

    bool project_3D(const libMesh::DenseMatrix<libMesh::Real> &polygon_1,
                    const libMesh::DenseMatrix<libMesh::Real> &polygon_2,
                    libMesh::DenseMatrix<libMesh::Real> &projection_1,
                    libMesh::DenseMatrix<libMesh::Real> &projection_2) {
        using namespace libMesh;

        typedef Intersector::Scalar Scalar;

        const int dim = 3;

        Scalar A[3 * 3], b[3];
        Scalar Ainv[3 * 3], binv[3];

        //      Scalar normal_master[3];
        //      Scalar normal_slave [3];

        Scalar isect_1[MAX_N_ISECT_POINTS * 3];
        Scalar isect_2[MAX_N_ISECT_POINTS * 3];

        libMesh::DenseMatrix<libMesh::Real> ref_polygon_1(polygon_1.m(), polygon_1.n());
        libMesh::DenseMatrix<libMesh::Real> ref_polygon_2(polygon_2.m(), polygon_2.n());
        libMesh::DenseMatrix<libMesh::Real> clipper_2(polygon_1.m(), 2);

        ref_polygon_1.resize(polygon_1.m(), polygon_1.n());
        clipper_2.resize(polygon_1.m(), dim - 1);
        ref_polygon_2.resize(polygon_2.m(), polygon_2.n());

        Intersector::triangle_make_affine_transform_3(&polygon_2.get_values()[0], A, b);
        Intersector::make_inverse_affine_transform_3(A, b, Ainv, binv);

        Intersector::apply_affine_transform_3(Ainv,
                                              binv,

                                              polygon_1.m(),
                                              &polygon_1.get_values()[0],
                                              &ref_polygon_1.get_values()[0]);

        // FIXME could store reference instead of computing it each time
        Intersector::apply_affine_transform_3(
            Ainv, binv, polygon_2.m(), &polygon_2.get_values()[0], &ref_polygon_2.get_values()[0]);

        for (uint i = 0; i < polygon_1.m(); ++i) {
            clipper_2(i, 0) = ref_polygon_2(i, 0);
            clipper_2(i, 1) = ref_polygon_2(i, 1);
        }

        const int n_projection_points = Intersector::project_surface_poly_onto_ref_poly(dim,
                                                                                        ref_polygon_1.m(),
                                                                                        &ref_polygon_1.get_values()[0],
                                                                                        polygon_2.m(),
                                                                                        &clipper_2.get_values()[0],
                                                                                        isect_1,
                                                                                        isect_2);

        // plot_polygon(clipper_2.n(),     clipper_2.m(),      &clipper_2.get_values()[0], "clipper");
        // plot_polygon(ref_polygon_1.n(), ref_polygon_1.m(),  &ref_polygon_1.get_values()[0], "clipped");
        // plot_polygon(3, n_projection_points,  isect_1, "clipper_isect");
        // plot_polygon(3, n_projection_points,  isect_2, "clipped_isect");

        if (!n_projection_points) {
            return false;
        }

        projection_1.resize(n_projection_points, 3);
        projection_2.resize(n_projection_points, 3);

        Intersector::apply_affine_transform_3(A, b, n_projection_points, isect_1, &projection_1.get_values()[0]);
        Intersector::apply_affine_transform_3(A, b, n_projection_points, isect_2, &projection_2.get_values()[0]);
        return true;
    }

    bool biorthgonal_weights(const int type, libMesh::Real &w_ii, libMesh::Real &w_ij) {
        using namespace libMesh;

        switch (type) {
            case EDGE2: {
                w_ii = 2.0;
                w_ij = -1.0;
                return true;
            }

            case TRI3: {
                w_ii = 3.0;
                w_ij = -1.0;
                return true;
            }

            case TET4: {
                w_ii = 4.0;
                w_ij = -1.0;
                return true;
            }
                // These do not work:
                // case QUAD4:
                // {
                //  w_ii = 4.0;
                //  w_ij = -1.0;
                //  return true;
                // }

                // case HEX8:
                // {
                //  w_ii = 8.0;
                //  w_ij = -1.0;
                //  return true;
                // }

            default: {
                w_ii = 1.0;
                w_ij = 0.0;

                static bool error_msg_printed = false;

                if (!error_msg_printed) {
                    std::cerr << "[Error] biorthgonal weights not supported for element type: " << type << std::endl;
                    error_msg_printed = true;
                }

                assert(false && "TODO: add the weights for the missing element");
                return false;
            }
        }
    }

    template <typename T>
    void convert_point_to_vector(const int dim, const libMesh::Point &point, std::vector<T> &point_vec) {
        point_vec.resize(dim);
        for (int i = 0; i < dim; ++i) {
            point_vec[i] = point(i);
        }
    }

    void mortar_normal_and_gap_assemble_weighted_biorth(const libMesh::FEBase &test_fe,
                                                        const int dim,
                                                        const libMesh::Point &surf_normal,
                                                        const libMesh::Point &plane_normal,
                                                        const libMesh::Real &plane_offset,
                                                        const libMesh::DenseMatrix<libMesh::Real> &weights,
                                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                                        libMesh::DenseVector<libMesh::Real> &gap,
                                                        const bool visdbg) {
        using namespace libMesh;

        if (normals.m() != test_fe.get_phi().size() || dim != normals.n()) {
            normals.resize(test_fe.get_phi().size(), dim);
            normals.zero();
            gap.resize(test_fe.get_phi().size());
            gap.zero();
        }

        const auto &test = test_fe.get_phi();
        // const auto &grad   = test_fe.get_dphi();
        const auto &point = test_fe.get_xyz();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_qp = test[0].size();

        std::vector<Real> surf_normal_v, plane_normal_v;
        convert_point_to_vector(dim, surf_normal, surf_normal_v);
        convert_point_to_vector(dim, plane_normal, plane_normal_v);

        DenseVector<Real> p(dim);
        // DenseMatrix<Real> v(dim); //visdbg

        for (uint i = 0; i < n_test; ++i) {
            for (uint qp = 0; qp < n_qp; ++qp) {
                p(0) = point[qp](0);
                p(1) = point[qp](1);

                if (dim > 2) {
                    p(2) = point[qp](2);
                }

                Real isect = 0;
                Intersector::intersect_ray_with_plane(
                    dim, 1, &p.get_values()[0], &surf_normal_v[0], &plane_normal_v[0], plane_offset, &isect);
                // assert(isect > 0);
                // if(visdbg) {
                //  v.get_values() = surf_normal_v; //visdbg
                //  v *= isect; //visdbg
                //  quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]); //visdbg
                // }

                auto biorth_test = weights(i, 0) * test[0][qp];

                for (uint k = 1; k < n_test; ++k) {
                    biorth_test += weights(i, k) * test[k][qp];
                }

                gap(i) += biorth_test * isect * JxW[qp];

                for (uint d = 0; d < dim; ++d) {
                    normals(i, d) += biorth_test * surf_normal(d) * JxW[qp];
                }
            }
        }
    }

    void mortar_normal_and_gap_assemble_weighted_biorth(const libMesh::FEVectorBase &test_fe,
                                                        const int dim,
                                                        const libMesh::Point &surf_normal,
                                                        const libMesh::Point &plane_normal,
                                                        const libMesh::Real &plane_offset,
                                                        const libMesh::DenseMatrix<libMesh::Real> &weights,
                                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                                        libMesh::DenseVector<libMesh::Real> &gap) {
        using namespace libMesh;

        if (normals.m() != test_fe.get_phi().size() / dim || dim != normals.n()) {
            normals.resize(test_fe.get_phi().size() / dim, dim);
            gap.resize(test_fe.get_phi().size());
        }

        normals.zero();
        gap.zero();

        const auto &test = test_fe.get_phi();
        const auto &point = test_fe.get_xyz();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_qp = test[0].size();

        DenseVector<Real> p(dim);
        DenseVector<Real> v(dim);

        DenseVector<Real> s_n(dim);
        DenseVector<Real> p_n(dim);

        for (uint i = 0; i < dim; ++i) {
            p_n(i) = plane_normal(i);
            s_n(i) = surf_normal(i);
        }

        for (uint qp = 0; qp < n_qp; ++qp) {
            p(0) = point[qp](0);
            p(1) = point[qp](1);

            if (dim > 2) {
                p(2) = point[qp](2);
            }

            Real isect = 0;

            Intersector::intersect_ray_with_plane(
                dim, 1, &p.get_values()[0], &s_n.get_values()[0], &p_n.get_values()[0], plane_offset, &isect);

            v = s_n;
            v *= isect;
            // quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]);

            for (uint i = 0; i < n_test; ++i) {
                auto biorth_test = weights(i, 0) * test[0][qp];

                for (uint k = 0; k < test.size(); ++k) {
                    biorth_test += weights(i, k) * test[k][qp];
                }

                gap(i) += biorth_test(0) * isect * JxW[qp];

                for (uint d = 0; d < dim; ++d) {
                    normals.get_values()[i] += biorth_test(d) * surf_normal(d) * JxW[qp];
                }
            }
        }

        // gap.print(std::cout);
    }

    void mortar_normal_and_gap_assemble(const uint dim,
                                        const libMesh::FEBase &test_fe,
                                        const libMesh::Point &surf_normal,
                                        const libMesh::Point &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap) {
        using namespace libMesh;
        DenseVector<Real> surf_normal_v(dim), plane_normal_v(dim);

        for (uint i = 0; i < dim; ++i) {
            surf_normal_v(i) = surf_normal(i);
            plane_normal_v(i) = plane_normal(i);
        }

        mortar_normal_and_gap_assemble(test_fe, surf_normal_v, plane_normal_v, plane_offset, normals, gap);
    }

    void mortar_normal_and_gap_assemble(const libMesh::FEBase &test_fe,
                                        const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                        const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap) {
        using namespace libMesh;

        const uint dim = plane_normal.size();

        if (normals.m() != test_fe.get_phi().size() || dim != normals.n()) {
            normals.resize(test_fe.get_phi().size(), dim);
            normals.zero();
            gap.resize(test_fe.get_phi().size());
            gap.zero();
        }

        const auto &test = test_fe.get_phi();
        // const auto &grad   = test_fe.get_dphi();
        const auto &point = test_fe.get_xyz();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_qp = test[0].size();

        DenseVector<Real> p(dim);

        for (uint i = 0; i < n_test; ++i) {
            for (uint qp = 0; qp < n_qp; ++qp) {
                p(0) = point[qp](0);
                p(1) = point[qp](1);

                if (dim > 2) {
                    p(2) = point[qp](2);
                }

                Real isect = 0;

                Intersector::intersect_ray_with_plane(dim,
                                                      1,
                                                      &p.get_values()[0],
                                                      &surf_normal.get_values()[0],
                                                      &plane_normal.get_values()[0],
                                                      plane_offset,
                                                      &isect);

                // printf("g: %g (%g, %g)\n", isect, p(1), plane_offset);
                // assert(isect > 0);

                gap(i) += test[i][qp] * isect * JxW[qp];

                for (uint d = 0; d < dim; ++d) {
                    normals(i, d) += test[i][qp] * surf_normal(d) * JxW[qp];
                }
            }
        }
    }

    void mortar_normal_and_gap_assemble(const uint dim,
                                        const libMesh::FEVectorBase &test_fe,
                                        const libMesh::Point &surf_normal,
                                        const libMesh::Point &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap) {
        using namespace libMesh;
        DenseVector<Real> surf_normal_v(dim), plane_normal_v(dim);

        for (uint i = 0; i < dim; ++i) {
            surf_normal_v(i) = surf_normal(i);
            plane_normal_v(i) = plane_normal(i);
        }

        mortar_normal_and_gap_assemble(test_fe, surf_normal_v, plane_normal_v, plane_offset, normals, gap);
    }

    void mortar_normal_and_gap_assemble(const libMesh::FEVectorBase &test_fe,
                                        const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                        const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap) {
        using namespace libMesh;

        const uint dim = plane_normal.size();

        if (normals.m() != test_fe.get_phi().size() / dim || dim != normals.n()) {
            normals.resize(test_fe.get_phi().size() / dim, dim);
            normals.zero();
            gap.resize(test_fe.get_phi().size());
            gap.zero();
        }

        const auto &test = test_fe.get_phi();
        const auto &point = test_fe.get_xyz();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_qp = test[0].size();

        DenseVector<Real> p(dim);
        DenseVector<Real> v(dim);

        for (uint qp = 0; qp < n_qp; ++qp) {
            p(0) = point[qp](0);
            p(1) = point[qp](1);

            if (dim > 2) {
                p(2) = point[qp](2);
            }

            Real isect = 0;

            Intersector::intersect_ray_with_plane(dim,
                                                  1,
                                                  &p.get_values()[0],
                                                  &surf_normal.get_values()[0],
                                                  &plane_normal.get_values()[0],
                                                  plane_offset,
                                                  &isect);

            v = surf_normal;
            v *= isect;
            // quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]);

            for (uint i = 0; i < n_test; ++i) {
                gap(i) += test[i][qp](0) * isect * JxW[qp];

                for (uint d = 0; d < dim; ++d) {
                    normals.get_values()[i] += test[i][qp](d) * surf_normal(d) * JxW[qp];
                }
            }
        }
    }

    // void mortar_normal_and_gap_assemble_biorth(const int type,
    //                                            const uint dim,
    //                                            const libMesh::FEVectorBase &test_fe,
    //                                            const libMesh::Point &surf_normal,
    //                                            const libMesh::Point &plane_normal,
    //                                            const libMesh::Real &plane_offset,
    //                                            const libMesh::DenseVector<libMesh::Real> &indicator,
    //                                            libMesh::DenseMatrix<libMesh::Real> &normals,
    //                                            libMesh::DenseVector<libMesh::Real> &gap) {
    //     using namespace libMesh;
    //     DenseVector<Real> surf_normal_v(dim), plane_normal_v(dim);

    //     for (uint i = 0; i < dim; ++i) {
    //         surf_normal_v(i) = surf_normal(i);
    //         plane_normal_v(i) = plane_normal(i);
    //     }

    //     mortar_normal_and_gap_assemble_biorth(
    //         type, test_fe, surf_normal_v, plane_normal_v, plane_offset, indicator, normals, gap);
    // }

    // bool intersect_2D(const libMesh::DenseMatrix<libMesh::Real> &poly1,
    //                   const libMesh::DenseMatrix<libMesh::Real> &poly2,
    //                   libMesh::DenseMatrix<libMesh::Real> &intersection) {
    //     double result_buffer[MAX_N_ISECT_POINTS * 2];
    //     int n_vertices_result;

    //     assert(Intersector::polygon_area_2(poly1.m(), &poly1.get_values()[0]) > 0);
    //     assert(Intersector::polygon_area_2(poly2.m(), &poly2.get_values()[0]) > 0);

    //     if (!
    //         // isector.
    //         intersect_convex_polygons(poly1.m(),
    //                                   &poly1.get_values()[0],
    //                                   poly2.m(),
    //                                   &poly2.get_values()[0],
    //                                   &n_vertices_result,
    //                                   result_buffer,
    //                                   DEFAULT_TOLLERANCE)) {
    //         return false;
    //     }

    //     assert(Intersector::polygon_area_2(n_vertices_result, result_buffer) > 0);

    //     intersection.resize(n_vertices_result, 2);
    //     std::copy(result_buffer, result_buffer + n_vertices_result * 2, &intersection.get_values()[0]);

    //     // plot_polygon(2, n_vertices_result, &intersection.get_values()[0], "isect/poly");
    //     return true;
    // }

    void make_polygon_from_quad4(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
        polygon.resize(e.n_nodes(), 2);

        for (int i = 0; i < e.n_nodes(); ++i) {
            for (int j = 0; j < 2; ++j) {
                polygon(i, j) = e.point(i)(j);
            }
        }
    }

    void make_polygon_from_high_order_quad(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
        polygon.resize(4, 2);

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 2; ++j) {
                polygon(i, j) = e.point(i)(j);
                // std::cout<<" polygon("<<i<<","<<j<<") = "<< polygon(i, j) <<std::endl;
            }
        }
    }

    void make_polygon_from_tri3(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
        polygon.resize(e.n_nodes(), 2);

        for (int i = 0; i < e.n_nodes(); ++i) {
            for (int j = 0; j < 2; ++j) {
                polygon(i, j) = e.point(i)(j);
            }
        }
    }

    template <int N>
    void discretize_segmented_curve(const double x[N],
                                    const double y[N],
                                    const int order,
                                    const libMesh::Elem &e,
                                    libMesh::DenseMatrix<libMesh::Real> &polygon) {
        using namespace libMesh;

        static const int Dim = 2;

        auto f = [&e](const double *x, double *fx) -> void {
            Point p;

            for (int d = 0; d < Dim; ++d) {
                p(d) = x[d];
            }

            Point fp = FE<Dim, libMesh::LAGRANGE>::map(&e, p);

            for (int d = 0; d < Dim; ++d) {
                fx[d] = fp(d);
            }
        };

        std::vector<double> all;
        std::vector<double> params_points, polyline;

        for (int i = 0; i < N; ++i) {
            const int ip1 = (i + 1) % N;
            const double from[Dim] = {x[i], y[i]};
            const double to[Dim] = {x[ip1], y[ip1]};

            discretize_curve<Dim>(f, from, to, order, params_points, polyline, 1e-6);
            all.insert(all.end(), polyline.begin(), polyline.end() - 2);
        }

        polygon.resize(all.size() / Dim, Dim);

        for (int i = 0; i < all.size(); i += Dim) {
            polygon(i / Dim, 0) = all[i];
            polygon(i / Dim, 1) = all[i + 1];
        }
    }

    void make_polygon_from_curved_tri6(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
        const double x[6] = {0, 0.5, 1, 0.5, 0, 0.0};
        const double y[6] = {0, 0.0, 0, 0.5, 1, 0.5};
        discretize_segmented_curve<6>(x, y, 2, e, polygon);
    }

    void make_polygon_from_curved_quad8(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
        using namespace libMesh;
        const double x[8] = {-1, 0, 1, 1, 1, 0, -1, -1};
        const double y[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
        discretize_segmented_curve<8>(x, y, 2, e, polygon);
    }

    void make_polygon_from_tri6(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
        polygon.resize(e.n_nodes() / 2, 2);
        for (int i = 0; i < e.n_nodes() / 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                polygon(i, j) = e.point(i)(j);
                // std::cout<<" polygon("<<i<<","<<j<<") = "<< polygon(i, j) <<std::endl;
            }
        }
    }

    void make_polyline(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polyline) {
        polyline.resize(2, 2);

        polyline(0, 0) = e.node_ref(0)(0);
        polyline(0, 1) = e.node_ref(0)(1);

        polyline(1, 0) = e.node_ref(1)(0);
        polyline(1, 1) = e.node_ref(1)(1);
    }

    void make_polygon(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
        // FIXME use libMesh enum types
        switch (e.n_nodes()) {
            case 2: {
                // works for lines too
                make_polygon_from_tri3(e, polygon);
                break;
            }
            case 3: {
                make_polygon_from_tri3(e, polygon);
                break;
            }

            case 4: {
                make_polygon_from_quad4(e, polygon);
                break;
            }

            case 6: {
                if (e.has_affine_map()) {
                    make_polygon_from_tri6(e, polygon);
                    // make_polygon_from_tri3(e, polygon);
                } else {
                    make_polygon_from_curved_tri6(e, polygon);
                }
                break;
            }

            case 8: {
                if (e.has_affine_map()) {
                    make_polygon_from_high_order_quad(e, polygon);
                    // make_polygon_from_quad4(e, polygon);
                } else {
                    make_polygon_from_curved_quad8(e, polygon);
                }
                break;
            }

            case 9: {
                // if(e.has_affine_map()) {
                make_polygon_from_high_order_quad(e, polygon);
                // make_polygon_from_quad4(e, polygon);
                // } else {
                //  make_polygon_from_curved_quad8(e, polygon);
                // }
                break;
            }

            default: {
                assert(false);
                break;
            }
        }
    }

    void make_polygon_3(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
        auto n_nodes = e.n_nodes();

        if (e.has_affine_map()) {
            if (is_tri(e.type())) {
                n_nodes = 3;
            } else if (is_quad(e.type())) {
                n_nodes = 4;
            } else {
                assert(false && "handle special case");
            }
        }

        polygon.resize(n_nodes, 3);

        for (int i = 0; i < n_nodes; ++i) {
            for (int j = 0; j < 3; ++j) {
                polygon(i, j) = e.point(i)(j);
            }
        }
    }

}  // namespace utopia

// clean-up
#ifdef USE_CLIPPER
#undef USE_CLIPPER
#endif

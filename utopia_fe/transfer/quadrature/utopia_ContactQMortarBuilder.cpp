#include "utopia_ContactQMortarBuilder.hpp"
#include "utopia_Intersect.hpp"
#include "utopia_SurfUtils.hpp"
#include "utopia_LibMeshShape.hpp"

namespace utopia {
    AffineContactQMortarBuilder3::AffineContactQMortarBuilder3(const Real search_radius)
    : trial_ir(3), test_ir(3), trial_box(3), test_box(3), total_intersection_volume(0.), search_radius(search_radius)
    {}

    bool AffineContactQMortarBuilder3::build(
            const Elem &trial,
            FEType trial_type,
            const int trial_side_num,
            const Elem &test,
            FEType test_type,
            const int test_side_num,
            QMortar &q_trial,
            QMortar &q_test)
    {
        auto trial_side = trial.build_side_ptr(trial_side_num);
        auto test_side  = test.build_side_ptr(test_side_num);

        compute_side_normal(3, *trial_side, trial_n);
        compute_side_normal(3, *test_side, test_n);

        const Real cos_angle = trial_n.contract(test_n);

        //if the angle is more than 60 degrees ( cos(60/180*pi) == 0.5 ) or has same orientation skip
        if(cos_angle >= -0.5) {
            return false;
        }

        trial_box.reset();
        test_box.reset();

        enlarge_box_from_side(3, *trial_side, trial_box, search_radius);
        enlarge_box_from_side(3, *test_side,  test_box,  search_radius);

        if(!trial_box.intersects(test_box)) return false;

        make_polygon_3(*trial_side, trial_polygon);
        make_polygon_3(*test_side,  test_polygon);

        if(!project_3D(trial_polygon,
                       test_polygon,
                       trial_isect,
                       test_isect)) {
            return false;
        }

        const Real area_slave = Intersector::polygon_area_3(test_polygon.m(), &test_polygon.get_values()[0]);
        const Real area   	  = Intersector::polygon_area_3(test_isect.m(),   &test_isect.get_values()[0]);

        // const Real relative_area = area/area_slave;
        const Real weight        = 1./area_slave;

        assert(area_slave > 0);
        assert(area > 0);
        assert(weight > 0);

        total_intersection_volume += area;

        const int order = order_for_l2_integral(3, trial, trial_type.order, test, test_type.order);

        make_composite_quadrature_on_surf_3D(trial_isect, weight, order, trial_ir);
        make_composite_quadrature_on_surf_3D(test_isect,  weight, order, test_ir);

        auto trial_trafo = std::make_shared<SideAffineTransform3>(trial, trial_side_num);
        auto test_trafo  = std::make_shared<SideAffineTransform3>(test, test_side_num);

        transform_to_reference_surf(*trial_trafo, trial.type(), trial_ir, q_trial);
        transform_to_reference_surf(*test_trafo,  test.type(),  test_ir,  q_test);

        const libMesh::Point pp = trial_side->point(0);
        const Real plane_offset = trial_n.contract(pp);


        auto n_qps = test_ir.n_points();
        gap_.resize(n_qps);
        normal_.resize(n_qps);

        double arr_trial_n[3] = { trial_n(0), trial_n(1), trial_n(2) };
        double arr_test_n[3]  = { test_n(0), test_n(1), test_n(2) };

        for(std::size_t qp = 0; qp < n_qps; ++qp) {
            Real isect = 0;
            const auto &p_temp = test_ir.get_points()[qp];
            double p[3] = { p_temp(0), p_temp(1), p_temp(2) };

            Intersector::intersect_ray_with_plane(
                3,
                1,
                p,
                arr_test_n,
                arr_trial_n,
                plane_offset,
                &isect);

            gap_[qp] = isect;
            normal_[qp] = test_n;
        }

        return true;
    }

    double AffineContactQMortarBuilder3::get_total_intersection_volume() const
    {
        return total_intersection_volume;
    }

    ///////////////////////////////////////////////////////////////////

    WarpedContactQMortarBuilder3::WarpedContactQMortarBuilder3(const Real search_radius)
    : trial_ir(3), test_ir(3), trial_box(3), test_box(3), total_intersection_volume(0.), search_radius(search_radius), ref_quad_order_(-1)
    {}

    bool WarpedContactQMortarBuilder3::build(
        const Elem &trial,
        FEType trial_type,
        const int trial_side_num,
        const Elem &test,
        FEType test_type,
        const int test_side_num,
        QMortar &q_trial,
        QMortar &q_test)
    {
        auto trial_side = trial.build_side_ptr(trial_side_num);
        auto test_side  = test.build_side_ptr(test_side_num);

        //compute avg normals
        SurfUtils::avg_normal(trial, trial_type, trial_n);
        SurfUtils::avg_normal(test,  test_type,  test_n);

        const Real cos_angle = trial_n.contract(test_n);

        //check angle
        //if the angle is more than 60 degrees ( cos(60/180*pi) == 0.5 ) or has same orientation skip
        if(cos_angle >= -0.5) {
            return false;
        }

        trial_box.reset();
        test_box.reset();

        enlarge_box_from_side(3, *trial_side, trial_box, search_radius);
        enlarge_box_from_side(3, *test_side,  test_box,  search_radius);

        if(!trial_box.intersects(test_box)) return false;

        //create polygons from warped surf . Options: 1) only inter nodes 2) discretized poly
        // using option 1

        const int order = order_for_l2_integral(3, trial, trial_type.order, test, test_type.order);
        init_ref_quad(order);

        make(*trial_side, trial_poly3_);
        make(*test_side, test_poly3_);

        const auto &p = test_side->node_ref(0);

        Plane3 plane = {
            {p(0), p(1), p(2)},
            {test_n(0), test_n(1), test_n(2)}
        };

        auto test_area = test_poly3_.area();

        //project on slave avg plane
        if(!project_intersect_and_map_quadrature(
            trial_poly3_,
            test_poly3_,
            plane,
            tri_q_points_,
            tri_q_weights_,
            1./test_area,
            composite_q_points_,
            composite_q_weights_
        )) {
            return false;
        }

        const bool use_newton = true;
        LibMeshShape<double, 3> trial_shape(*trial_side, trial_type, use_newton);
        LibMeshShape<double, 3> test_shape(*test_side, test_type, use_newton);

        //create quad points
        bool ok = trial_shape.make_quadrature(
            plane.n,
            composite_q_points_,
            composite_q_weights_,
            q_trial,
            trial_gap_
        ); assert(ok);

        if(!ok) return false;

        ok = test_shape.make_quadrature(
            plane.n,
            composite_q_points_,
            composite_q_weights_,
            q_test,
            test_gap_
        ); assert(ok);

        if(ok) {
            auto n_qps = test_gap_.size();
            gap_.resize(n_qps);
            normal_.resize(n_qps);

            for(std::size_t qp = 0; qp < n_qps; ++qp) {
                gap_[qp] = trial_gap_[qp] + test_gap_[qp];
                normal_[qp] = test_n;
            }
        }

        return ok;
    }

    double WarpedContactQMortarBuilder3::get_total_intersection_volume() const
    {
        return total_intersection_volume;
    }

    void WarpedContactQMortarBuilder3::init_ref_quad(const int order)
    {
        if(ref_quad_order_ != order) {
            libMesh::QGauss ir(2, libMesh::Order(order));
            ir.init(libMesh::TRI3);

            auto n_qps = ir.n_points();
            tri_q_points_.resize(n_qps);
            tri_q_weights_.resize(n_qps);

            for(int k = 0; k < n_qps; ++k) {
                auto &qp = ir.get_points()[k];

                for(int d = 0; d < 3; ++d) {
                    tri_q_points_[k][d] = qp(d);
                }

                tri_q_weights_[k] = ir.get_weights()[k];
            }
        }
    }

}

#include "utopia_ContactQMortarBuilder.hpp"
#include "utopia_Intersect.hpp"
#include "utopia_SurfUtils.hpp"


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

		const Real relative_area = area/area_slave;
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
		return true;
	}

	double AffineContactQMortarBuilder3::get_total_intersection_volume() const
	{
		return total_intersection_volume;
	}

	///////////////////////////////////////////////////////////////////

	WarpedContactQMortarBuilder3::WarpedContactQMortarBuilder3(const Real search_radius)
	: trial_ir(3), test_ir(3), trial_box(3), test_box(3), total_intersection_volume(0.), search_radius(search_radius)
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
		make_polygon_3(*trial_side, trial_polygon);
		make_polygon_3(*test_side,  test_polygon);

		const auto &p = test_side->node_ref(0);

		Plane3 plane = {
			{p(0), p(1), p(2)},
			{test_n(0), test_n(1), test_n(2)}
		};

		// std::vector<Polygon3::Vector> composite_q_points;
		// std::vector<Polygon3::Scalar> composite_q_weights;

		// bool ok = project_intersect_and_map_quadrature(
		// 	trial_polygon,
		// 	test_polygon,
		// 	plane,
		// 	//ref-quad-rule
		// 	{{1./3., 1./3.}},
		// 	{1.},
		// 	1.,
		// 	composite_q_points,
		// 	composite_q_weights
		// ); utopia_test_assert(ok);



		// ok = left_shape.make_quadrature(
		//     plane.n,
		//     composite_q_points,
		//     composite_q_weights,
		//     left_q
		// ); utopia_test_assert(ok);

		// ok = right_shape.make_quadrature(
		//     plane.n,
		//     composite_q_points,
		//     composite_q_weights,
		//     right_q
		// ); utopia_test_assert(ok);


		//project on slave avg plane
		//create quad points
		return false;
	}

	double WarpedContactQMortarBuilder3::get_total_intersection_volume() const
	{
		return total_intersection_volume;
	}

}

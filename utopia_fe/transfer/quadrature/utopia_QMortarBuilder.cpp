#include "utopia_QMortarBuilder.hpp"
#include "MortarAssemble.hpp"

namespace utopia {
	bool QMortarBuilder2::build(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type,
		QMortar &q_trial,
		QMortar &q_test)
	{
		make_polygon(trial, trial_pts);
		make_polygon(test, test_pts);

		if(intersect_2D(trial_pts, test_pts, intersection)) {
			total_intersection_volume += fabs(isector.polygon_area_2(intersection.m(), &intersection.get_values()[0]));
			const libMesh::Real weight = isector.polygon_area_2(test_pts.m(), &test_pts.get_values()[0]);

			const int order = order_for_l2_integral(2, trial, trial_type.order, test, test_type.order);
			make_composite_quadrature_2D(intersection, weight, order, composite_ir);
			auto trial_trans  = std::make_shared<AffineTransform2>(trial);
			auto test_trans   = std::make_shared<AffineTransform2>(test);

			transform_to_reference(*trial_trans,  trial.type(), composite_ir, q_trial);


			// if(vol2surf) {
			// 	transform_to_reference_surf(*dest_trans, dest_el.type(), composite_ir, dest_ir);
			// } else {
			transform_to_reference(*test_trans,   test.type(),  composite_ir, q_test);
			// }

			return true;
		} else {
			return false;
		}
	}

	bool QMortarBuilder3::build(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type,
		QMortar &q_trial,
		QMortar &q_test)
	{
		make_polyhedron(trial, trial_poly);
		make_polyhedron(test,  test_poly);

		if(intersect_3D(trial_poly, test_poly, intersection)) {
			total_intersection_volume += compute_volume(intersection);
			const libMesh::Real weight = compute_volume(test_poly);
			auto trial_trans = std::make_shared<AffineTransform3>(trial);
			std::shared_ptr<Transform> test_trans;

			bool vol2surf = false;

			if(is_tri(test.type()) || is_quad(test.type())) {
				test_trans = std::make_shared<Transform2>(test);
				vol2surf = true;
			} else {
				test_trans = std::make_shared<AffineTransform3>(test);
			}

			const int order = order_for_l2_integral(2, trial, trial_type.order, test, test_type.order);

			if(vol2surf) {
				shell_poly.resize(intersection.n_nodes, 3);
				std::copy(intersection.points, intersection.points + intersection.n_nodes * intersection.n_dims, &shell_poly.get_values()[0]);
				make_composite_quadrature_on_surf_3D(shell_poly, 1./weight, order, composite_ir);

			// plot_polygon(3, test_poly.n_nodes, test_poly.points, "polygon/" + std::to_string(comm.rank()) + "/p");
			} else {
				make_composite_quadrature_3D(intersection, weight, order, composite_ir);
			}

			transform_to_reference(*trial_trans, trial.type(), composite_ir,  q_trial);

			if(vol2surf) {
				transform_to_reference_surf(*test_trans, test.type(), composite_ir, q_test);
			} else {
				transform_to_reference(*test_trans, test.type(), composite_ir,  q_test);
			}

			return true;
		} else {
			return false;
		}
	}
}
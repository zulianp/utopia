#include "utopia_QMortarBuilder.hpp"
#include "MortarAssemble.hpp"

#include <cmath>

namespace utopia {
	using Scalar   = double;
	// using SizeType = int;
	using Vector2  = Intersector::Vector2;


	static const int INSIDE  = 1;
	static const int ON_EDGE = 2;
	static const int OUTSIDE = 0;

	inline bool intersect_convex_polygon_with_polyline(
		const SizeType n_vertices_1,
		const Scalar *polygon_1,
		const SizeType n_vertices_2,
		const Scalar *polygon_2,
		SizeType *n_vertices_result,
		Scalar *result,
		const Scalar tol)
	{
		short is_inside    [MAX_N_ISECT_POINTS];
		Scalar input_buffer[MAX_N_ISECT_POINTS * 2];

		const SizeType n_elements_2 = n_vertices_2 * 2;

		SizeType n_input_points  = n_vertices_2;
		SizeType n_output_points = 0;


		Scalar * input  = input_buffer;
		Scalar * output = result;

		generic_copy(n_elements_2, polygon_2, input);


		Vector2 e1, e2, p, s, e;

		for(SizeType i = 0; i < n_vertices_1; ++i) {
			if(i > 0) {
				generic_ptr_swap(Scalar, input, output);
				n_input_points = n_output_points;
				n_output_points = 0;
			}

			const SizeType i2x 	= i * 2;
			const SizeType i2y  = i2x + 1;

			const SizeType i2p1x = 2 * (((i + 1) == n_vertices_1)? 0 : (i + 1));
			const SizeType i2p1y = i2p1x + 1;

			e1 = Intersector::vec_2(polygon_1[i2x], polygon_1[i2y]);
			e2 = Intersector::vec_2(polygon_1[i2p1x], polygon_1[i2p1y]);

			SizeType n_outside = 0;
			for(SizeType j = 0; j < n_input_points; ++j) {
				const SizeType jx = j * 2;
				const SizeType jy = jx + 1;

				p = Intersector::vec_2(input[jx], input[jy]);
				is_inside[j] = Intersector::inside_half_plane(e1, e2, p, tol);
				n_outside += is_inside[j] == OUTSIDE;
			}

			if(n_input_points - n_outside == 0) return false;

			for(SizeType j = 0; j < n_input_points; ++j) {
				const SizeType jx = j * 2;
				const SizeType jy = jx + 1;

				const SizeType jp1 = (j + 1 == n_input_points)? 0 : (j + 1);
				const SizeType jp1x = jp1 * 2;
				const SizeType jp1y = jp1x + 1;

				s = Intersector::vec_2(input[jx], input[jy]);
				e = Intersector::vec_2(input[jp1x], input[jp1y]);

				if(is_inside[j]) {
					n_output_points = Intersector::append_point(n_output_points, s, output);

					if( ( is_inside[j] != ON_EDGE ) && ( !is_inside[jp1] ) ) {
						n_output_points = Intersector::append_point(n_output_points, Intersector::intersect_lines(e1, e2, s, e), output);
					}
				} else if(is_inside[jp1]) {
					n_output_points = Intersector::append_point(n_output_points, Intersector::intersect_lines(e1, e2, s, e), output);
				} else if(is_inside[j] == ON_EDGE && is_inside[jp1] == ON_EDGE) {
					n_output_points = Intersector::append_point(n_output_points, s, output);
					n_output_points = Intersector::append_point(n_output_points, e, output);
				}
			}

			if(n_output_points < 2) {
				return false;
			}

			n_output_points = Intersector::collapse_quasi_equal_points(n_output_points, output, tol);

			if(n_output_points < 2) {
				return false;
			}
		}

		*n_vertices_result = n_output_points;

		if(output != result) {
			const SizeType n_intersection_elements = n_output_points * 2;
			generic_copy(n_intersection_elements, output, result);
		}

		return n_output_points >= 2;
	}

	inline Scalar polyline_length(const SizeType n_vertices, const Scalar * points)
	{
		Scalar res = 0.;
		for(SizeType i = 1; i < n_vertices; ++i) {
			const SizeType ix_prev = (i-1) * 2;
			const SizeType iy_prev = (i-1) * 2 + 1;

			const SizeType ix = i * 2;
			const SizeType iy = i * 2 + 1;

			const Scalar dx = points[ix] - points[ix_prev];
			const Scalar dy = points[iy] - points[iy_prev];

			res += std::sqrt(dx*dx + dy*dy);
		}

		return res;
	}


	bool QMortarBuilder2::build(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type,
		QMortar &q_trial,
		QMortar &q_test)
	{
		if(test.dim() == 1) {
			//volume to surface

			make_polygon(trial, trial_pts);
			make_polyline(test, test_pts);

			std::vector<double> intersection_temp(test_pts.get_values().size()*2, 0.);


			SizeType n_res_pts;
			if(intersect_convex_polygon_with_polyline(
					trial_pts.get_values().size()/2,
					&trial_pts.get_values()[0],
					test_pts.get_values().size()/2,
					&test_pts.get_values()[0],
					&n_res_pts,
					&intersection_temp[0],
					1e-10)) {

				intersection_temp.resize(n_res_pts * 2);
				intersection.resize(n_res_pts, 2);
				intersection.get_values() = intersection_temp;

				total_intersection_volume += polyline_length(n_res_pts, &intersection.get_values()[0]);
				const Scalar weight = polyline_length(test_pts.n(), &test_pts.get_values()[0]);

				const int order = order_for_l2_integral(2, trial, trial_type.order, test, test_type.order);

				auto trial_trans  = std::make_shared<AffineTransform2>(trial);
				auto test_trans   = std::make_shared<Transform1>(test);

				make_composite_quadrature_on_surf_2D(intersection, 1./weight, order, composite_ir);

				transform_to_reference(*trial_trans, trial.type(), composite_ir, q_trial);
				transform_to_reference_surf(*test_trans, test.type(), composite_ir, q_test);
				return true;
			}

			return false;
		}


		make_polygon(trial, trial_pts);
		make_polygon(test, test_pts);

		if(intersect_2D(trial_pts, test_pts, intersection)) {
			total_intersection_volume += fabs(Intersector::polygon_area_2(intersection.m(), &intersection.get_values()[0]));
			const libMesh::Real weight = Intersector::polygon_area_2(test_pts.m(), &test_pts.get_values()[0]);

			const int order = order_for_l2_integral(2, trial, trial_type.order, test, test_type.order);
			make_composite_quadrature_2D(intersection, weight, order, composite_ir);
			auto trial_trans  = std::make_shared<AffineTransform2>(trial);
			auto test_trans   = std::make_shared<AffineTransform2>(test);

			transform_to_reference(*trial_trans,  trial.type(), composite_ir, q_trial);
			transform_to_reference(*test_trans, test.type(), composite_ir, q_test);

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

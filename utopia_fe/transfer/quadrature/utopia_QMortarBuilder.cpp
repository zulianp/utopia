#include "utopia_QMortarBuilder.hpp"
#include "MortarAssemble.hpp"
#include "utopia_Socket.hpp"
#include "utopia_libmesh_Utils.hpp"

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


	bool QMortarBuilder1::build(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type,
		QMortar &q_trial,
		QMortar &q_test)
	{
		assert(trial.n_nodes() == 2);
		assert(test.n_nodes()  == 2);

		const auto &p1 = trial.node_ref(0);
		const auto &p2 = trial.node_ref(1);

		const auto &q1 = test.node_ref(0);
		const auto &q2 = test.node_ref(1);

		u = p2 - p1;
		v = q2 - q1;

		auto len_u = u.norm();
		auto len_v = v.norm();

		u /= len_u;
		v /= len_v;

		const double cos_angle = u * v;

		if(std::abs(std::abs(cos_angle) - 1.) > 1e-14) {
			//not collinear
			return false;
		}

		w = q2 - p1;
		auto len_w = w.norm();

		if(len_w != 0.) {
			w /= len_w;
			const double cos_angle_2 = w * u;
			if(std::abs(std::abs(cos_angle_2) - 1.) > 1e-14) {
				//not on the same plane
				return false;
			}
		}

		for(int i = 0; i < LIBMESH_DIM; ++i) {
			min_p(i) = std::min(p1(i), p2(i));
			max_p(i) = std::max(p1(i), p2(i));
			
			min_q(i) = std::min(q1(i), q2(i));
			max_q(i) = std::max(q1(i), q2(i));

			intersection[0](i) = std::max(min_p(i), min_q(i));
			intersection[1](i) = std::min(max_p(i), max_q(i));
		}

		r = intersection[1] - intersection[0];

		auto isect_len = r.norm_sq();
		if(isect_len < 1e-16) {
			return false;
		}

		libMesh::DenseMatrix<libMesh::Real> line(2, LIBMESH_DIM);

		for(int i = 0; i < LIBMESH_DIM; ++i) {
			line(0, i) = intersection[0](i);
			line(1, i) = intersection[1](i);
		}

		const int order = order_for_l2_integral(1, trial, trial_type.order, test, test_type.order);
		make_composite_quadrature_on_surf_2D(line, 1./len_v, order, composite_ir);
		total_intersection_volume += std::sqrt(isect_len);

		auto trial_trans = std::make_shared<Transform1>(trial);
		auto test_trans  = std::make_shared<Transform1>(test);

		transform_to_reference(*trial_trans, trial.type(), composite_ir, q_trial);
		transform_to_reference(*test_trans,  test.type(),  composite_ir, q_test);
		return true;
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
				transform_to_reference(*test_trans,  test.type(),  composite_ir, q_test);
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

	static double counter_weight(libMesh::ElemType type)
	{
		if(is_tri(type)) {
			return 1./2.;
		} else if(is_quad(type)) {
			return 1./2.;
		} else {
			assert(false && "add special case");
			return 1.;
		}
	}

	bool QMortarBuilderShell2::build(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type,
		QMortar &q_trial,
		QMortar &q_test) 
	{
		compute_normal(trial, trial_normal);
		compute_normal(test,  test_normal);

		auto angle = trial_normal * test_normal;

		if(std::abs(angle - 1.) > 1e-14) {
			//not coplanar and not same orientation
			return false;
		}

		make_polygon_3(trial, trial_pts);
		make_polygon_3(test,  test_pts);

		if(!intersect()) {
			return false;
		}

		const int order = order_for_l2_integral(2, trial, trial_type.order, test, test_type.order);
		const Scalar weight = Intersector::polygon_area_3(test_pts.m(), &test_pts.get_values()[0]);
		
		make_composite_quadrature_on_surf_3D(intersection, 1./weight, order, composite_ir);

		total_intersection_volume += Intersector::polygon_area_3(intersection.m(), &intersection.get_values()[0]);
		
		auto trial_trans = std::make_shared<Transform2>(trial);
		auto test_trans  = std::make_shared<Transform2>(test);

		transform_to_reference(*trial_trans, trial.type(), composite_ir, q_trial);
		transform_to_reference(*test_trans,  test.type(),  composite_ir, q_test);

		double cw = counter_weight(test.type());
		for(auto &w : q_test.get_weights()) {
			w *= cw;
		}

		cw = counter_weight(trial.type());
		for(auto &w : q_trial.get_weights()) {
			w *= cw;
		}
		
		// assert(false);

		// plot_polygon(3, test_pts.m(), 	  &test_pts.get_values()[0], 	 "test");
		// plot_polygon(3, trial_pts.m(),    &trial_pts.get_values()[0], 	 "trial");

		// static int n_isect = 0;
		// plot_polygon(3, intersection.m(), &intersection.get_values()[0], "isect" + std::to_string(n_isect) + "/poly");
		// plot_quad_points(3, composite_ir.get_points(), "isect" + std::to_string(n_isect) + "/qp");
		// plot_quad_points(3, q_trial.get_points(), "qptrial");
		// ++n_isect;

		return true;
	}

	bool QMortarBuilderShell2::intersect()
	{
		using namespace libMesh;
		typedef Intersector::Scalar Scalar;

		ref_trial_pts.resize(trial_pts.m(), trial_pts.n());
		ref_trial_pts_2.resize(trial_pts.m(), 2);

		ref_test_pts.resize(test_pts.m(), test_pts.n());
		ref_test_pts_2.resize(test_pts.m(), 2);


		Intersector::triangle_make_affine_transform_3(&test_pts.get_values()[0], A, b);
		Intersector::make_inverse_affine_transform_3(A, b, Ainv, binv);

		Intersector::apply_affine_transform_3(Ainv, binv,
										 	  trial_pts.m(),
										      &trial_pts.get_values()[0],
										      &ref_trial_pts.get_values()[0]
										     );

		Intersector::apply_affine_transform_3(Ainv, binv, test_pts.m(),
										 	  &test_pts.get_values()[0],
										      &ref_test_pts.get_values()[0]
										      );

		for(uint i = 0; i < ref_trial_pts.m(); ++i) {
			ref_trial_pts_2(i, 0) = ref_trial_pts(i, 0);
			ref_trial_pts_2(i, 1) = ref_trial_pts(i, 1);
			assert(approxeq(0., ref_trial_pts(i, 2)));
		}

		for(uint i = 0; i < ref_test_pts.m(); ++i) {
			ref_test_pts_2(i, 0) = ref_test_pts(i, 0);
			ref_test_pts_2(i, 1) = ref_test_pts(i, 1);
			assert(approxeq(0., ref_test_pts(i, 2)));
		}

		if(!intersect_2D(ref_trial_pts_2, ref_test_pts_2, ref_intersection_2)) {
			return false;
		}

		ref_intersection_slave.resize(ref_intersection_2.m(),  3);
		ref_intersection_master.resize(ref_intersection_2.m(), 3);
		intersection.resize(ref_intersection_2.m(), 3);

		for(uint i = 0; i < ref_intersection_2.m(); ++i) {
			ref_intersection_slave(i, 0) = ref_intersection_2(i, 0);
			ref_intersection_slave(i, 1) = ref_intersection_2(i, 1);
			ref_intersection_slave(i, 2) = 0.;
		}

		Intersector::apply_affine_transform_3(
			A, b,
			ref_intersection_slave.m(),
			&ref_intersection_slave.get_values()[0],
			&intersection.get_values()[0]
		);

		return true;
	}

	class Plane3 {
	public:
		using Vector = Intersector::Vector3;
		using Scalar = Intersector::Scalar;

		inline Scalar signed_dist(const Vector &q) const
		{
			return dot(n, q - p);
		}

		Vector n;
		Vector p;
	};

	class HalfSpace3 {
	public:
		using Vector = Intersector::Vector3;
		using Scalar = Intersector::Scalar;

		Vector n;
		Scalar d;

		//negative is inside the half space
		inline Scalar signed_dist(const Vector &q) const
		{
			return dot(n, q) - d;
		}
	};

	class Line3 {
	public:
		using Vector = Intersector::Vector3;
		using Scalar = Intersector::Scalar;

		Vector p0;
		Vector p1;
	};

	class Polygon3 {
	public:
		using Vector = Intersector::Vector3;
		using Scalar = Intersector::Scalar;

		inline void clear()
		{
			points.clear();
		}

		std::vector<Vector> points;
	};

	class HPolyhedron3 {
	public:
		std::vector<HalfSpace3> half_spaces;
	};

	static int intersect_planes(const Plane3 &plane_1, const Plane3 &plane_2, Line3 &L)
	{
	    using Vector = Plane3::Vector;
	    using Scalar = Plane3::Scalar;
	    
	    Vector   u  = cross(plane_1.n, plane_2.n);          // cross product
	    Scalar   ax = (u.x >= 0 ? u.x : -u.x);
	    Scalar   ay = (u.y >= 0 ? u.y : -u.y);
	    Scalar   az = (u.z >= 0 ? u.z : -u.z);

	    // test if the two planes are parallel
	    if ((ax+ay+az) < 1e-16) {        // plane_1 and plane_2 are near parallel
	        // test if disjoint or coincide
	        Vector v = plane_2.p -  plane_1.p;
	        if (dot(plane_1.n, v) == 0)          // plane_2.p lies in plane_1
	            return 1;                    // plane_1 and plane_2 coincide
	        else 
	            return 0;                    // plane_1 and plane_2 are disjoint
	    }

	    // plane_1 and plane_2 intersect in a line
	    // first determine max abs coordinate of cross product
	    int      maxc;                       // max coordinate
	    if (ax > ay) {
	        if (ax > az)
	             maxc =  1;
	        else maxc = 3;
	    }
	    else {
	        if (ay > az)
	             maxc =  2;
	        else maxc = 3;
	    }

	    // next, to get a point on the intersect line
	    // zero the max coord, and solve for the other two
	    Vector  iP;                // intersect point
	    Scalar    d1, d2;            // the constants in the 2 plane equations
	    d1 = -dot(plane_1.n, plane_1.p);  // note: could be pre-stored  with plane
	    d2 = -dot(plane_2.n, plane_2.p);  // ditto

	    switch (maxc) {             // select max coordinate
	    case 1:                     // intersect with x=0
	        iP.x = 0;
	        iP.y = (d2*plane_1.n.z - d1*plane_2.n.z) /  u.x;
	        iP.z = (d1*plane_2.n.y - d2*plane_1.n.y) /  u.x;
	        break;
	    case 2:                     // intersect with y=0
	        iP.x = (d1*plane_2.n.z - d2*plane_1.n.z) /  u.y;
	        iP.y = 0;
	        iP.z = (d2*plane_1.n.x - d1*plane_2.n.x) /  u.y;
	        break;
	    case 3:                     // intersect with z=0
	        iP.x = (d2*plane_1.n.y - d1*plane_2.n.y) /  u.z;
	        iP.y = (d1*plane_2.n.x - d2*plane_1.n.x) /  u.z;
	        iP.z = 0;
	    }

	    L.p0 = iP;
	    L.p1 = iP + u;
	    return 2;
	}

	static int intersect(const Line3 &line, const Plane3 &plane, Line3::Scalar &t, const Line3::Scalar tol = 1e-10)
	{
		using Vector = Line3::Vector;
		using Scalar = Line3::Scalar;

		auto ray_dir = line.p1 - line.p0;
		auto len = length(ray_dir);
		ray_dir /= len;

		const Scalar cos_angle = -dot(plane.n, ray_dir);

		if(fabs(cos_angle) < tol) {
			t = 0;
			return 0;
		}

		//distance from plane
		const Scalar dist = plane.signed_dist(line.p1);

		//point in ray trajectory
		t = dist/cos_angle;

		if(t >= -tol && t <= 1. + tol) {
			return 1;
		} else {
			return -1;
		}
	}

	static int intersect(const Line3 &line, const HalfSpace3 &half_space, Line3 &result, const Line3::Scalar tol = 1e-10)
	{
		using Vector = Line3::Vector;
		using Scalar = Line3::Scalar;

		auto d0 = half_space.signed_dist(line.p0);
		auto d1 = half_space.signed_dist(line.p1);

		if(std::signbit(d0) && std::signbit(d1)) {
			//inside half-space 
			result = line;
			return 2;
		}

		if(!std::signbit(d0) && !std::signbit(d1)) {
			if(d0 > 0 || d1 > 0) {
				//no intersection or point intersection
				return 0;
			}
		}

		if(approxeq(d0, 0., tol) && approxeq(d1, 0., tol)) {
			//aligned with half-space boundary
			result = line;
			return 3;
		}

		Plane3 plane;
		plane.n = half_space.n;
		plane.p = half_space.n * half_space.d;

		Scalar t = 0.;
		auto code = intersect(line, plane, t, tol);

		if(d0 < 0) {
			result.p0 = line.p0;
			result.p1 = line.p0 + t * (line.p1 - line.p0);
		} else if(d1 < 0) {
			result.p0 = line.p0 + t * (line.p1 - line.p0);
			result.p1 = line.p1;
		}

		return 1;
	}

	static int intersect(const Polygon3 &polygon, const HPolyhedron3 &h_poly, Polygon3 &result, const Scalar tol = 1e-10)
	{
		using Vector = Polygon3::Vector;
		using Scalar = Polygon3::Scalar;

		const std::size_t n_points 	  = polygon.points.size();
		const std::size_t n_half_spaces = h_poly.half_spaces.size();

		// std::vector<Scalar> signed_dist(n_half_spaces * n_points, 0.);

		// std::size_t n_all_negative = 0;
		// for(std::size_t i = 0; i < n_half_spaces; ++i) {
		// 	auto offset_i = i * n_points;
			
		// 	bool all_negative = true;
		// 	bool all_positive = true;
			
		// 	for(std::size_t j = 0; j < n_points; ++j) {
		// 		auto d = h_poly.half_spaces[i].signed_dist(poly.points[j]);
		// 		all_positive = all_positive && d > 0.;
		// 		all_negative = all_negative && d < 0.;
		// 		signed_dist[offset_i + j] = d;
		// 	}

		// 	if(all_positive) {
		// 		//points are all outside 
		// 		return 0;
		// 	}

		// 	n_all_negative += all_negative;
		// }

		// if(n_all_negative == n_half_spaces) {
		// 	//points are completely contained in convex-set
		// 	result = poly;
		// 	return 1;
		// }

		// static const short LINE_OUTSIDE     = 0;
		// static const short LINE_INSIDE      = 1;
		// static const short LINE_COINCIDE 	= 2;
		// static const short LINE_INTERSECT_1 = 3;
		// static const short LINE_INTERSECT_2 = 4;

		// std::vector<short> intersection_types(n_points, LINE_OUTSIDE);
		// std::vector<std::vector<std::size_t>> intersecting_with(n_points);
		// std::vector<short> temp(n_half_spaces, LINE_OUTSIDE);

		// // //for each line find all intersections
		// for(std::size_t i = 0; i < n_points; ++i) {
		// 	const auto ip1 = (i + 1 == n_points) ? 0 : (i+1);

		// 	for(std::size_t k = 0; k < n_half_spaces; ++k) {
		// 		const auto offset_k = k * n_points;
				
		// 		const auto d1 = signed_dist[offset_k + i];
		// 		const auto d2 = signed_dist[offset_k + i];

		// 		if(d1 < 0 && d2 < 0) {
		// 			temp[k] = LINE_INSIDE;
		// 		} else if(d1 > 0 && d2 > 0) {
		// 			temp[k] = LINE_OUTSIDE;
		// 		} else if(std::signbit(d1) != std::signbit(d2)) {
		// 			temp[k] = LINE_INTERSECT_1;
		// 		} else {
		// 			assert(approxeq(d1, 0.) && approxeq(d2, 0.));
		// 			temp[k] = LINE_COINCIDE;
		// 		}
		// 	}

		// 	for(std::size_t k = 0; k < n_half_spaces; ++k) {
		// 		if(temp[k] == LINE_INTERSECT_1) {
		// 			intersecting_with[i].push_back(k);

		// 			if(intersection_types[i] == LINE_INTERSECT_1) {
		// 				intersection_types[i] = LINE_INTERSECT_2;
		// 				continue;
		// 			}
		// 		}

		// 		intersection_types[i] = std::max(temp[k], intersection_types[i]);
		// 	}
		// }

		result.clear();

		for(std::size_t i = 0; i < n_points; ++i) {
			const auto ip1 = (i + 1 == n_points) ? 0 : (i+1);
			Line3 line;
			line.p0 = polygon.points[i];
			line.p1 = polygon.points[ip1];

			int n_intersected = 0;
			bool line_outside_h_poly = false;
			bool aligned_with_separating_plane = false;
			Line3 isect_line;
			
			for(std::size_t k = 0; k < n_half_spaces; ++k) {
				auto ret = intersect(line, h_poly.half_spaces[k], isect_line, tol);
				switch(ret) {
					case 0: {
						//outside half-space
						line_outside_h_poly = true;
						break;
					}

					case 1: {
						//intersection with separating plane
						line = isect_line;
						++n_intersected;
						break;
					}

					case 2: {
						//completely inside
						// line3 = isect_line; //they are the same
						break;
					}

					case 3: {
						aligned_with_separating_plane = true;
						break;
					}
				}

				if(line_outside_h_poly) { break; }
			}

			if(line_outside_h_poly) {
				//the segment is to be skipped
				continue;
			}

			if(n_intersected) {

			}
		}

		return true;
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

		static int n_isect = 0;

		if(n_isect == 20) {
			std::cout << "HERE";
		}

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

				// if(n_isect == 791) {
					plot_polygon(3, shell_poly.m(), &shell_poly.get_values()[0], "isect" + std::to_string(n_isect) + "/poly");
				// }


				n_isect++;
			// plot_polygon(3, test_poly.n_nodes, test_poly.points, "polygon/" + std::to_string(comm.rank()) + "/p");
			} else {
				make_composite_quadrature_3D(intersection, weight, order, composite_ir);
			}

			transform_to_reference(*trial_trans, trial.type(), composite_ir,  q_trial);

			if(vol2surf) {
				// transform_to_reference_surf(*test_trans, test.type(), composite_ir, q_test);
				transform_to_reference(*test_trans, test.type(), composite_ir,  q_test);

				double cw = counter_weight(test.type());
				for(auto &w : q_test.get_weights()) {
					w *= cw;
				}



			} else {
				transform_to_reference(*test_trans, test.type(), composite_ir,  q_test);
			}

			return true;
		} else {
			return false;
		}
	}
}

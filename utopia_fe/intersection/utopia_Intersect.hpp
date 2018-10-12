#ifndef UTOPIA_INTERSECT_HPP
#define UTOPIA_INTERSECT_HPP

#include "utopia_intersector.hpp"
#include "utopia_Socket.hpp"

#include <libmesh/elem.h>

namespace utopia {

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

		inline Vector barycenter() const
		{
			Vector res = {0., 0., 0.};

			for(const auto &p : points) {
				res += p;
			}

			res /= points.size();
			return res;
		}

		const Vector &point(const SizeType i) const
		{
			return points[i];
		}

		inline Scalar area() const
		{
			auto n = points.size();

			if(n == 3) {
				return 0.5 * Intersector::trapezoid_area_3(points[0], points[1], points[2]);
			}

			Vector b = barycenter();

			Scalar result = 0;
			for(SizeType i = 0; i < n; ++i) {
				const SizeType ip1 = (i + 1 == n)? 0 : (i + 1);
				result += Intersector::trapezoid_area_3(b, points[i], points[ip1]);
			}

			return 0.5 * result;
		}

		inline void remove_duplicate_points(const Scalar tol) {
			if(points.empty()) return;

			const SizeType n_original = points.size();
			const SizeType last = n_original - 1;
			std::vector<bool> keep(n_original, true);

			SizeType n = n_original;
			
			for(SizeType i = 1; i < n_original; ++i) {
				const SizeType i_prev = i - 1;
					
				auto d = distance(points[i_prev], points[i]);

				if(d <= tol) {
					n--;
					keep[i] = false;
				} else {
					keep[i] = true;
				}
			}

			if(keep[last]) {
				auto d = distance(points[0], points[last]);
				if(d <= tol) {
					n--;
					keep[last] = false;
				} else {
					keep[last] = true;
				}
			}

			if(n == points.size()) return;

			std::vector<Vector> old_points = std::move(points);
			points.reserve(n);

			for(std::size_t i = 0; i < n_original; ++i) {
				if(keep[i]) {
					points.push_back(old_points[i]);
				}
			}
		}

		inline void resize(const std::size_t new_size)
		{
			points.resize(new_size);
		}

		inline void plot(const std::string &name) const {

			std::vector<double> xyz;
			xyz.reserve(points.size() * 3);

			for(const auto &p : points) {
				xyz.push_back(p.x);
				xyz.push_back(p.y);
				xyz.push_back(p.z);
			}

			plot_polygon(3, points.size(), &xyz[0], name);
		}

		inline bool equals(const Polygon3 &other, const Scalar tol = 1e-16) const
		{
			std::size_t n = points.size();
			if(n != other.points.size()) return false;

			for(std::size_t i = 0; i < n; ++i) {
				auto d = length(points[i] - other.points[i]);
				if(d > tol) { return false; }
			}

			return true;
		}

		std::vector<Vector> points;
	};

	class HPolyhedron3 {
	public:
		using Vector = HalfSpace3::Vector;
		using Scalar = HalfSpace3::Scalar;

		std::vector<HalfSpace3> half_spaces;
	};

	inline int intersect_planes(const Plane3 &plane_1, const Plane3 &plane_2, Line3 &L)
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

	inline int intersect(const Line3 &line, const Plane3 &plane, Line3::Scalar &t, const Line3::Scalar tol = 1e-10)
	{
		using Vector = Line3::Vector;
		using Scalar = Line3::Scalar;

		auto ray_dir = line.p1 - line.p0;
		auto len = length(ray_dir);
		ray_dir /= len;

		const Scalar cos_angle = dot(plane.n, ray_dir);

		if(fabs(cos_angle) < tol) {
			t = 0;
			return 0;
		}

		//distance from plane
		const Scalar dist = plane.signed_dist(line.p0);

		//point in ray trajectory
		t = dist/cos_angle;
		t /= len;

		if(dist <= 0. && cos_angle > 0.) {
			//intersection with ray coming from behind the plane
			t *= -1.;
		} else if(dist > 0. && cos_angle < 0.) {
			//intersection with ray coming from the front of the plane
			t *= -1;
		}

		if(t >= -tol && t <= 1. + tol) {
			//inside segment
			return 1;
		} else {
			//outside segment
			return -1;
		}
	}

	inline int intersect(const Line3 &line, const HalfSpace3 &half_space, Line3 &result, const Line3::Scalar tol = 1e-10)
	{
		using Vector = Line3::Vector;
		using Scalar = Line3::Scalar;

		auto d0 = half_space.signed_dist(line.p0);
		auto d1 = half_space.signed_dist(line.p1);

		if(std::signbit(d0) && std::signbit(d1)) {
			//inside half-space 
			result = line;
			return 3;
		}

		if(!std::signbit(d0) && !std::signbit(d1)) {
			if(d0 > tol || d1 > tol) {
				//no intersection or point intersection
				return 0;
			}
		}

		if(approxeq(d0, 0., tol) && approxeq(d1, 0., tol)) {
			//aligned with half-space boundary
			result = line;
			return 4;
		}

		Plane3 plane;
		plane.n = half_space.n;
		plane.p = half_space.n * half_space.d;

		Scalar t = 0.;
		auto code = intersect(line, plane, t, tol);

		if(code == 0) {
			//line is collinear
			assert(false);
		}

		//no actual intersection with segment
		if(code == -1) {
			if(d0 <= tol || d1 <= tol) {
				result = line;
				return 3;
			} else {
				return 0;
			}
		}

		assert(t >= -tol);
		assert(t <= 1. + tol);

		auto u = line.p1 - line.p0;

		if(d0 < d1) {
			result.p0 = line.p0;
			result.p1 = line.p0 + t * u;
			return 1;
		} else {
			result.p0 = line.p0 + t * u;
			result.p1 = line.p1;
			return 2;
		}
	}

	inline int intersect(const Polygon3 &polygon, const HPolyhedron3 &h_poly, Polygon3 &result, const Polygon3::Scalar tol = 1e-10)
	{
		using Vector = Polygon3::Vector;
		using Scalar = Polygon3::Scalar;

		Line3 line, isect_line;

		const std::size_t n_half_spaces = h_poly.half_spaces.size();
		Polygon3 in = polygon;

		for(std::size_t k = 0; k < n_half_spaces; ++k) {
			const std::size_t n_points = in.points.size();
			result.clear();

			for(std::size_t i = 0; i < n_points; ++i) {
				const auto ip1 = (i + 1 == n_points) ? 0 : (i+1);

				line.p0 = in.points[i];
				line.p1 = in.points[ip1];

				int n_intersected = 0;
				bool line_outside_h_poly = false;
				bool aligned_with_separating_plane = false;

				bool orginal_vertex[2] = { true, true };

				auto ret = intersect(line, h_poly.half_spaces[k], isect_line, tol);
				
				switch(ret) {
					case 0: {
						//outside half-space
						line_outside_h_poly = true;
						break;
					}
					case 1: 
					{
						orginal_vertex[1] = false;
						line = isect_line;
						++n_intersected;
						break;
					}
					case 2:
					{
						//intersection with separating plane
						orginal_vertex[0] = false;
						line = isect_line;
						++n_intersected;
						break;
					}

					case 3: {
						//completely inside
						// line3 = isect_line; //they are the same
						break;
					}

					case 4: {
						aligned_with_separating_plane = true;
						break;
					}
				}

				
				if(line_outside_h_poly) {
					//the segment is to be skipped and not to be used
					continue;
				}

				result.points.push_back(isect_line.p0);

				if(n_intersected && !orginal_vertex[1]) {
					result.points.push_back(isect_line.p1);
				}
			}

			result.remove_duplicate_points(tol);
			if(result.points.size() < 3) return 0;

			in = result;
		}

		return (result.points.size() >= 3)? 1 : 0;
	}

	inline void make(const libMesh::Point &p, Polygon3::Vector &q)
	{
		q.x = p(0);
		q.y = p(1);
		q.z = p(2);
	}

	inline void make(const libMesh::Elem &elem, Polygon3 &polygon) {
		assert(is_tri(elem.type()) || is_quad(elem.type()));
		auto n_nodes = is_tri(elem.type()) ? 3 : 4;

		polygon.resize(n_nodes);
		
		for(auto i = 0; i < n_nodes; ++i) {
			const auto &p = elem.node_ref(i);
			auto &q = polygon.points[i];
			make(p, q);
		}
	}

	inline void make(
		const HPolyhedron3::Vector &p0,
		const HPolyhedron3::Vector &p1,
		const HPolyhedron3::Vector &p2,
		HalfSpace3 &h
		)
	{
		h.n = cross(p1 - p0, p2 - p0);
		h.n /= length(h.n);
		h.d = dot(h.n, p0);
	}

	inline void make(const libMesh::Elem &elem, HPolyhedron3 &h_poly)
	{
		using Vector = HPolyhedron3::Vector;

		Vector p0, p1, p2, p3;
		Vector u, v;

		make(elem.node_ref(0), p0);
		make(elem.node_ref(1), p1);
		make(elem.node_ref(2), p2);
		make(elem.node_ref(3), p3);

		if(is_tet(elem.type())) {
			h_poly.half_spaces.resize(4);

			//face 0 = [0, 1, 3]
			make(p0, p1, p3, h_poly.half_spaces[0]);

			//face 1 = [1, 2, 3]
			make(p1, p2, p3, h_poly.half_spaces[1]);

			//face 2 = [0, 3, 2]
			make(p0, p3, p2, h_poly.half_spaces[2]);

			//face 3 = [1, 2, 0]
			make(p1, p2, p0, h_poly.half_spaces[3]);

			//fix orientation
			Vector b = p0; b += p1; b += p2; b += p3;
			b /= 4.0;

			for(int i = 0; i < 4; ++i) {
				if(h_poly.half_spaces[i].signed_dist(b) > 0.) {
					//change face orientation
					h_poly.half_spaces[i].n *= -1.;
					h_poly.half_spaces[i].d *= -1.;
				}
			}
		
		} else if(is_hex(elem.type())) {
			h_poly.half_spaces.resize(6);

			Vector p4, p5, p6, p7;

			make(elem.node_ref(4), p4);
			make(elem.node_ref(5), p5);
			make(elem.node_ref(6), p6);
			make(elem.node_ref(7), p7);

			//face 0 = [0, 1, 5, 4]
			make(p0, p1, p5, h_poly.half_spaces[0]);

			//face 1 = [1, 2, 6, 5]
			make(p1, p2, p6, h_poly.half_spaces[1]);

			//face 2 = [3, 7, 6, 2]
			make(p3, p7, p6, h_poly.half_spaces[2]);

			//face 3 = [0, 4, 7, 3]
			make(p0, p4, p7, h_poly.half_spaces[3]);

			//face 4 = [2, 1, 0, 3]
			make(p2, p1, p0, h_poly.half_spaces[4]);

			//face 5 = [6, 7, 4, 5]
			make(p6, p7, p4, h_poly.half_spaces[5]);

		} else {
			assert(false);
		}
	}

}

#endif //UTOPIA_INTERSECT_HPP

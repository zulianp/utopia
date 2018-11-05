#include "utopia_Project.hpp"
#include "utopia_HouseholderTransformation.hpp"
#include "utopia_triangulate.hpp"

namespace utopia {
	template<class Point3, class Point2>
	inline static void project_points(
		const HouseholderTransformation &trafo,
		const Point3 &translation,
		const std::vector<Point3> &in,
		std::vector<Point2> &out)
	{
		auto n_points = in.size();
		out.resize(n_points);

		Point3 p_in, p_out;
		for(std::size_t i = 0; i < n_points; ++i) {

			//translate to origin
			p_in = in[i] - translation;
			trafo.apply(p_in, p_out);

			//remove z coord
			out[i].x = p_out.x;
			out[i].y = p_out.y;
			// assert(approxeq(p_out.z, 0.)); //only if coplanar
		}
	}

	template<class Point2, class Point3>
	inline static void unproject_points(
		const HouseholderTransformation &trafo,
		const Point3 &translation,
		const std::vector<Point2> &in,
		std::vector<Point3> &out)
	{
		Point3 in_p;

		auto n_points = in.size();
		out.resize(n_points);

		for(std::size_t i = 0; i < n_points; ++i) {
			const auto &q = in[i];
			in_p.x = q.x;
			in_p.y = q.y;
			in_p.z = 0.;
			trafo.apply(in_p, out[i]);

			out[i] += translation;
		}
	}

	bool project_and_intersect(
		const Polygon3 &poly1,
		const Polygon3 &poly2,
		const Plane3 &plane,
		std::vector<Polygon3> &result
	)
	{
		using Vector3 = Polygon3::Vector;
		using Vector2 = Polygon2::Vector;

		auto v = plane.n;
		v.z += 1.;
		auto len = length(v);
		v /= len;

		HouseholderTransformation trafo(v);

		Polygon2 projected_1, projected_2;
		project_points(trafo, plane.p, poly1.points, projected_1.points);
		project_points(trafo, plane.p, poly2.points, projected_2.points);

		const auto n_points_1 = poly1.points.size();
		const auto n_points_2 = poly2.points.size();

		std::vector<Polygon2> result_2d;
		if(!Polygon2::intersect(projected_1, projected_2, result_2d)) {
			return false;
		}

		const auto n_islands = result_2d.size();
		result.resize(n_islands);

		for(std::size_t k = 0; k < n_islands; ++k) {
			unproject_points(trafo, plane.p, result_2d[k].points, result[k].points);
		}

		return true;
	}

	bool project_intersect_and_map_quadrature(
		const Polygon3 &poly1,
		const Polygon3 &poly2,
		const Plane3 &plane,
		//ref-quad-rule
		const std::vector<Polygon2::Vector> &q_points,
		const std::vector<Polygon2::Scalar> &q_weights,
		const Polygon3::Scalar weight,
		std::vector<Polygon3::Vector>  &composite_q_points,
		std::vector<Polygon3::Scalar>  &composite_q_weights
	)
	{
		using Vector3 = Polygon3::Vector;
		using Vector2 = Polygon2::Vector;

		auto v = plane.n;
		v.z += 1.;
		auto len = length(v);
		v /= len;

		HouseholderTransformation trafo(v);

		Polygon2 projected_1, projected_2;
		project_points(trafo, plane.p, poly1.points, projected_1.points);
		project_points(trafo, plane.p, poly2.points, projected_2.points);

		const auto n_points_1 = poly1.points.size();
		const auto n_points_2 = poly2.points.size();

		std::vector<Polygon2> result_2d;
		if(!Polygon2::intersect(projected_1, projected_2, result_2d)) {
			return false;
		}

		const auto n_islands = result_2d.size();
		std::vector<std::vector<int>> triangulations(n_islands);

		for(std::size_t k = 0; k < n_islands; ++k) {
			triangulate_polygon(result_2d[k], triangulations[k]);
		}

		std::vector<Polygon2::Vector> composite_q_points_2d;
		
		map_quadrature_rule(
			q_points,
			q_weights,
			weight,
			result_2d,
			triangulations,
			composite_q_points_2d,
			composite_q_weights
		);

		auto n_qps = composite_q_points_2d.size();
		unproject_points(trafo, plane.p, composite_q_points_2d, composite_q_points);
		return true;
	}

	void map_quadrature_rule(
		const std::vector<Polygon2::Vector> &q_points,
		const std::vector<Polygon2::Scalar> &q_weights,
		const Polygon3::Scalar weight,
		const std::vector<Polygon2> &domain_of_integration,
		const std::vector<std::vector<int>> &triangulations,
		std::vector<Polygon2::Vector>  &composite_q_points,
		std::vector<Polygon2::Scalar>  &composite_q_weights
	)
	{
		using Vector2 = Polygon2::Vector;

		std::size_t n_qps_rule = q_points.size();
		std::size_t n_qps = 0;
		std::size_t n_islands = domain_of_integration.size();

		for(std::size_t k = 0; k < n_islands; ++k) {
			n_qps += triangulations[k].size() * n_qps_rule;
		}

		composite_q_points.resize(n_qps);
		composite_q_weights.resize(n_qps);

		Vector2 u, v;

		std::size_t qp = 0;
		for(std::size_t k = 0; k < n_islands; ++k) {
			const auto &tri = triangulations[k];
			const auto &points = domain_of_integration[k].points;

			const auto n_triangles = tri.size()/3;

			for(std::size_t i = 0; i < n_triangles; ++i) {
				const auto i3 = i * 3;
				const auto &o = points[tri[i3]];
				
				u = points[tri[i3 + 1]] - o;
				v = points[tri[i3 + 2]] - o;

				for(std::size_t j = 0; j < n_qps_rule; ++j) {
					composite_q_points[qp] = o;
					composite_q_points[qp] += q_points[j].x * u;
					composite_q_points[qp] += q_points[j].y * v;
					++qp;
				}
			}
		}
	}
}

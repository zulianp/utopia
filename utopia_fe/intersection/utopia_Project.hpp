#ifndef UTOPIA_PROJECT_HPP
#define UTOPIA_PROJECT_HPP

#include "utopia_Intersect.hpp"
#include "utopia_Polygon2.hpp"
#include <vector>

namespace utopia {

	// namespace new_types {

		template<typename Scalar, int Dim>
		using Vector = moonolith::Vector<Scalar, Dim>;

		// template<typename Scalar, int Dim>
		// class Transform {
		// public:
		// 	using Vector = utopia::new_types::Vector<Scalar, Dim>;

		// 	//FIXME use custom type instead of stl vector
		// 	using Matrix = std::vector<Scalar>;

		// 	bool apply(const Vector &ref, Vector &world) const = 0;
		// 	bool apply_inverse(const Vector &ref, Vector &world) const = 0;
		// 	bool jacobian(const Vector &ref, Matrix &mat) const = 0;
		// };

		template<typename Scalar, int Dim>
		class Ray {
		public:
			Vector<Scalar, Dim> o, dir;
		};

		template<typename Scalar, int Dim>
		class Shape {
		public:
			virtual ~Shape() {}

			// Transform<Scalar, Dim> &transform();
			virtual bool intersect(
				const Ray<Scalar, Dim> &ray,
				Scalar &t) = 0;
		};

		template<typename Scalar, int Dim>
		class LibMeshShape final {
		public:
			LibMeshShape(const libMesh::Elem &elem)
			: elem_(elem)
			{}

			bool intersect(
				const Ray<Scalar, Dim> &ray,
				Scalar &t) override
			{
				return false;
			}

		private:
			const libMesh::Elem &elem_;
		};

		

	// }

	/**
	 * @brief poly1 and poly2 can be non-planar and non-convex polygons
	 */
	bool project_and_intersect(
		const Polygon3 &poly1,
		const Polygon3 &poly2,
		const Plane3 &plane,
		std::vector<Polygon3> &result
	);

	void map_quadrature_rule(
		const std::vector<Polygon2::Vector> &q_points,
		const std::vector<Polygon2::Scalar> &q_weights,
		const Polygon3::Scalar weight,
		const std::vector<Polygon2> &domain_of_integration,
		const std::vector<std::vector<int>> &triangulations,
		std::vector<Polygon2::Vector>  &composite_q_points,
		std::vector<Polygon2::Scalar>  &composite_q_weights
	);

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
	);

}

#endif //UTOPIA_PROJECT_HPP

#ifndef UTOPIA_PROJECT_HPP
#define UTOPIA_PROJECT_HPP

#include "utopia_Intersect.hpp"
#include "utopia_Polygon2.hpp"


#include <vector>



namespace utopia {

        template<typename Scalar, int Dim>
        using Vector = moonolith::Vector<Scalar, Dim>;

        // template<typename Scalar, int Dim>
        // class Ray {
        // public:
        //     Vector<Scalar, Dim> o, dir;
        // };

        // template<typename Scalar, int Dim>
        // class Shape {
        // public:
        //     virtual ~Shape() {}

        //     virtual bool intersect(
        //         const Ray<Scalar, Dim> &ray,
        //         Scalar &t) = 0;
        // };

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

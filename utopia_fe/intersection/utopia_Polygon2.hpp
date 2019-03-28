#ifndef UTOPIA_POLYGON_2_HPP
#define UTOPIA_POLYGON_2_HPP

#include "utopia_Intersect.hpp"
#include <vector>

namespace utopia {

    class Polygon2 {
    public:
        using Vector = Intersector::Vector2;
        using Scalar = Intersector::Scalar;

        inline std::size_t size() const
        {
            return points.size();
        }

        inline bool empty() const
        {
            return points.empty();
        }

        Scalar area() const;

        static bool intersect(const Polygon2 &poly1, const Polygon2 &poly2, std::vector<Polygon2> &result, const Scalar tol = DEFAULT_TOLLERANCE);
        std::vector<Vector> points;
    };

}

#endif //UTOPIA_POLYGON_2_HPP

#include "utopia_Polygon2.hpp"
#include <cmath>
#include "clipper.hpp"

namespace utopia {

    bool Polygon2::intersect(const Polygon2 &poly1,
                             const Polygon2 &poly2,
                             std::vector<Polygon2> &result,
                             Polygon2::Scalar tol) {
        using namespace ClipperLib;

        auto n_vertices_1 = poly1.size();
        auto n_vertices_2 = poly2.size();

        double min_x = poly1.points[0].x;
        double min_y = poly1.points[0].y;
        double max_x = poly1.points[0].x;
        double max_y = poly1.points[0].y;

        for (int i = 1; i < n_vertices_1; ++i) {
            min_x = std::min(min_x, poly1.points[i].x);
            max_x = std::max(max_x, poly1.points[i].x);

            min_y = std::min(min_y, poly1.points[i].y);
            max_y = std::max(max_y, poly1.points[i].y);
        }

        for (int i = 1; i < n_vertices_2; ++i) {
            min_x = std::min(min_x, poly2.points[i].x);
            max_x = std::max(max_x, poly2.points[i].x);

            min_y = std::min(min_y, poly2.points[i].y);
            max_y = std::max(max_y, poly2.points[i].y);
        }

        const double denom = std::max(max_x - min_x, max_y - min_y);
        const double cut_off = std::min(1e16, 1e16 / denom);  // 1e18 is the max represented value

        Paths subj(1), clip(1), solution;

        subj[0].reserve(n_vertices_1);
        clip[0].reserve(n_vertices_2);

        for (int i = 0; i < n_vertices_1; ++i) {
            subj[0].push_back(IntPoint((poly1.points[i].x - min_x) * cut_off, (poly1.points[i].y - min_y) * cut_off));
        }

        for (int i = 0; i < n_vertices_2; ++i) {
            clip[0].push_back(IntPoint((poly2.points[i].x - min_x) * cut_off, (poly2.points[i].y - min_y) * cut_off));
        }

        // perform intersection ...
        Clipper c;
        c.AddPaths(subj, ptSubject, true);
        c.AddPaths(clip, ptClip, true);
        c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);

        if (solution.empty() || solution[0].size() < 3) return false;

        assert(solution.size() == 1);

        auto n_islands = solution.size();
        result.resize(n_islands);

        std::size_t n_actual_intersection = 0;

        for (std::size_t k = 0; k < n_islands; ++k) {
            auto n_nodes = solution[k].size();
            if (n_nodes < 3) {
                continue;
            }

            result[n_actual_intersection].points.resize(n_nodes);

            for (std::size_t i = 0; i < n_nodes; ++i) {
                result[n_actual_intersection].points[i].x = solution[k][i].X / cut_off + min_x;
                result[n_actual_intersection].points[i].y = solution[k][i].Y / cut_off + min_y;
            }

            ++n_actual_intersection;
        }

        if (n_actual_intersection == 0) {
            result.clear();
            return false;
        }

        result.resize(n_actual_intersection);

        return true;
    }

    Polygon2::Scalar Polygon2::area() const {
        const auto n = size();
        std::vector<Scalar> arr(n * 2);

        std::size_t index = 0;
        for (std::size_t i = 0; i < n; ++i, index += 2) {
            arr[index] = points[i].x;
            arr[index + 1] = points[i].y;
        }

        return Intersector::polygon_area_2(n, &arr[0]);
    }
}  // namespace utopia

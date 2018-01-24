
#include "clipper.hpp"
#include <iostream>

#include "utopia_intersector.hpp"

// #define USE_CLIPPER 1

namespace utopia {
	bool intersect_convex_polygons(const int n_vertices_1,
								   const double * polygon_1,
								   const int n_vertices_2,
								   const double * polygon_2,
								   int *n_vertices_result,
								   double *result_buffer,
								   double tol) {

		// static bool lib_msg = false;

#ifdef USE_CLIPPER
		// if(!lib_msg) {
		// 	std::cout << "using clipper" << std::endl;
		// 	lib_msg = true;
		// }

		using namespace ClipperLib;

		double min_x = polygon_1[0];
		double min_y = polygon_1[1];
		double max_x = polygon_1[0];
		double max_y = polygon_1[1];

		for(int i = 1; i < n_vertices_1; ++i) {
			const int i2 = i * 2;
			min_x = std::min(min_x, polygon_1[i2]);
			max_x = std::max(max_x, polygon_1[i2]);

			min_y = std::min(min_y, polygon_1[i2+1]);
			max_y = std::max(max_y, polygon_1[i2+1]);
		}

		for(int i = 0; i < n_vertices_2; ++i) {
			const int i2 = i * 2;
			min_x = std::min(min_x, polygon_2[i2]);
			max_x = std::max(max_x, polygon_2[i2]);

			min_y = std::min(min_y, polygon_2[i2+1]);
			max_y = std::max(max_y, polygon_2[i2+1]);
		}


		const double cut_off = 1e10/std::max(max_x - min_x, max_y - min_y); //1e18 is the max represented value 

		Paths subj(1), clip(1), solution;

		subj[0].reserve(n_vertices_1);
		clip[0].reserve(n_vertices_2);

		for(int i = 0; i < n_vertices_1; ++i) {
			const int i2 = i * 2;
			subj[0].push_back(IntPoint((polygon_1[i2] - min_x) * cut_off , (polygon_1[i2+1] - min_y) * cut_off));
		}

		for(int i = 0; i < n_vertices_2; ++i) {
			const int i2 = i * 2;
			clip[0].push_back(IntPoint((polygon_2[i2] - min_x) * cut_off, (polygon_2[i2+1] - min_y) * cut_off));
		}

		//perform intersection ...
		Clipper c;
		c.AddPaths(subj, ptSubject, true);
		c.AddPaths(clip, ptClip, true);
		c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);

		if(solution.empty() || solution[0].size() < 3) return false;

		assert(solution.size() == 1);

		*n_vertices_result = solution[0].size();

		for(std::size_t i = 0; i < solution[0].size(); ++i) {
			const int i2 = i * 2;
			result_buffer[i2] = solution[0][i].X/cut_off + min_x;
			result_buffer[i2 + 1] = solution[0][i].Y/cut_off + min_y;
		}
		
		return true;

#else //USE_CLIPPER

		// if(!lib_msg) {
		// 	std::cout << "using homemade" << std::endl;
		// 	lib_msg = true;
		// }


		Intersector isector;
		return isector.intersect_convex_polygons(n_vertices_1, polygon_1, n_vertices_2, polygon_2, n_vertices_result, result_buffer, tol);
#endif //USE_CLIPPER
	}
}

//clean-up
#ifdef USE_CLIPPER
#undef USE_CLIPPER
#endif

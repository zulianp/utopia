
#include "clipper.hpp"
#include <iostream>

#include "utopia_intersector.hpp"

// #define USE_CLIPPER

namespace utopia {
	bool intersect_convex_polygons(const int n_vertices_1,
		const double * polygon_1,
		const int n_vertices_2,
		const double * polygon_2,
		int *n_vertices_result,
		double *result_buffer,
		double tol) {

		static bool lib_msg = false;

#ifdef USE_CLIPPER
		if(!lib_msg) {
			std::cout << "using clipper" << std::endl;
			lib_msg = true;
		}

		using namespace ClipperLib;

		const double cut_off = 1e12;

		Paths subj(1), clip(1), solution;

		subj[0].reserve(n_vertices_1);
		clip[0].reserve(n_vertices_2);

		for(int i = 0; i < n_vertices_1; ++i) {
			const int i2 = i * 2;
			subj[0].push_back(IntPoint(polygon_1[i2] * cut_off, polygon_1[i2+1] * cut_off));
		}

		for(int i = 0; i < n_vertices_2; ++i) {
			const int i2 = i * 2;
			clip[0].push_back(IntPoint(polygon_2[i2] * cut_off, polygon_2[i2+1] * cut_off));
		}

	//perform intersection ...
		Clipper c;
		c.AddPaths(subj, ptSubject, true);
		c.AddPaths(clip, ptClip, true);
		c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);

		if(solution.empty()) return false;

		assert(solution.size() == 1);

		*n_vertices_result = solution[0].size();

		for(std::size_t i = 0; i < solution[0].size(); ++i) {
			const int i2 = i * 2;
			result_buffer[i2] = solution[0][i].X/cut_off;
			result_buffer[i2 + 1] = solution[0][i].Y/cut_off;
		}
		
		return solution[0].size() >= 3;

#else //USE_CLIPPER

		if(!lib_msg) {
			std::cout << "using homemade" << std::endl;
			lib_msg = true;
		}


		Intersector isector;
		return isector.intersect_convex_polygons(n_vertices_1, polygon_1, n_vertices_2, polygon_2, n_vertices_result, result_buffer, tol);
#endif //USE_CLIPPER
	}
}

//clean-up
#ifdef USE_CLIPPER
#undef USE_CLIPPER
#endif

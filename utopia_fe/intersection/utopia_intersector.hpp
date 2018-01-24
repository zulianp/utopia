#ifndef UTOPIA_INTERSECTOR_HPP
#define UTOPIA_INTERSECTOR_HPP

#include <math.h>
#include "opencl_adapter.hpp"

#define USE_DOUBLE_PRECISION
#define DEFAULT_TOLLERANCE 1e-12

namespace utopia {

	class Intersector : public moonolith::OpenCLAdapter {
	public:
		//I do not know why the compiler wants this...
		template<typename T>
		inline static T sqrt(const T v) {
			return std::sqrt(v);
		}

		#include "all_kernels.cl"
	};

	typedef utopia::Intersector::PMesh Polyhedron;
	typedef utopia::Intersector::FEObject FEObject;

	bool intersect_convex_polygons(const int n_vertices_1,
								   const double * polygon_1,
								   const int n_vertices_2,
								   const double * polygon_2,
								   int *n_vertices_result,
								   double *result_buffer,
								   double tol);
}

#undef mortar_assemble

#endif //UTOPIA_INTERSECTOR_HPP

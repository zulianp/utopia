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
}

#undef mortar_assemble

#endif //UTOPIA_INTERSECTOR_HPP

#ifndef UTOPIA_ELASTIC_MATERIAL_HPP
#define UTOPIA_ELASTIC_MATERIAL_HPP

#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_LameeParameters.hpp"

namespace utopia {

	template<class Matrix, class Vector>
	class ElasticMaterial {
	public:
		virtual ~ElasticMaterial() {}
		// virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) = 0;
		virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) = 0;

		virtual bool stress(const Vector &x, Vector &result) {
			assert(false && "implement me");
			return false;
		}

		virtual void clear() {}
	};

}

#endif //UTOPIA_ELASTIC_MATERIAL_HPP

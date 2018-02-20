#ifndef UTOPIA_FUNCTION_MATERIAL_HPP
#define UTOPIA_FUNCTION_MATERIAL_HPP 

#include "utopia_ElasticMaterial.hpp"
#include "utopia_Function.hpp"

namespace utopia {
	template<class Matrix, class Vector>
	class FunctionMaterial : public ElasticMaterial<Matrix, Vector> {
	public:
		FunctionMaterial(const std::shared_ptr< Function<Matrix, Vector> > &fun)
		: fun_(fun)
		{}

		virtual ~FunctionMaterial() {}

		virtual bool assemble_hessian_and_gradient(Vector &x, Matrix &hessian, Vector &gradient)
		{
			bool ok = fun_->gradient(x, gradient);
			ok 	   &= fun_->hessian(x, hessian);
			return ok;
		}
		
		virtual void clear() {}

	private:
		std::shared_ptr< Function<Matrix, Vector> > fun_;
	};
}

#endif //UTOPIA_FUNCTION_MATERIAL_HPP

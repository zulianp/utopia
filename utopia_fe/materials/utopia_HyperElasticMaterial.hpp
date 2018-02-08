#ifndef UTOPIA_HYPER_ELASTIC_MATERIAL_HPP
#define UTOPIA_HYPER_ELASTIC_MATERIAL_HPP

namespace utopia {
	template<class Matrix, class Vector>
	class HyperElasticMaterial {
	public:
		virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) = 0;
	};
}


#endif //UTOPIA_HYPER_ELASTIC_MATERIAL_HPP

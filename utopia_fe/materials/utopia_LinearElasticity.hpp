#ifndef UTOPIA_LINEAR_ELASTICITY_HPP
#define UTOPIA_LINEAR_ELASTICITY_HPP

namespace utopia {

	template<class FunctionSpaceT, class Matrix, class Vector>
	class LinearElasticity : public ElasticMaterial<Matrix, Vector> {
	public:
		LinearElasticity(FunctionSpaceT &V, const LameeParameters &params)
		: V_(V), params_(params), initialized_(false)
		{}

		bool init(Matrix &hessian)
		{
			if(initialized_) return true;

			auto u = trial(V_);
			auto v = test(V_);

			auto mu     = params_.var_mu();
			auto lambda = params_.var_lambda();

			auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) ); 
			auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

			auto b_form = integral((2. * mu) * inner(e_u, e_v) + lambda * inner(div(u), div(v)));
			
			initialized_ = assemble(b_form, hessian);
			return initialized_;
		}

		// bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
		bool assemble_hessian_and_gradient(Vector &x, Matrix &hessian, Vector &gradient) override
		{
			if(!init(hessian)) {
				return false;
			}

			gradient = hessian * x;
			return true;
		}

		void clear() override
		{
			initialized_ = false;
		}

	private:
		FunctionSpaceT &V_;
		LameeParameters params_;
		bool initialized_;
	};
}

#endif //UTOPIA_LINEAR_ELASTICITY_HPP

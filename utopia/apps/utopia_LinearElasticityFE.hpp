#ifndef UTOPIA_LINEAR_ELASTICITY_FE_HPP
#define UTOPIA_LINEAR_ELASTICITY_FE_HPP

#include "utopia_Core.hpp"
#include "utopia_PetscDM.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_petsc.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_petsc_dma_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_PhaseField.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_Input.hpp"


namespace utopia {

	template<class FunctionSpace>
	class LinearElasticityFE final : public Operator<typename FunctionSpace::Vector>, public Configurable {
	public:
	    // using Mesh             = typename FunctionSpace::Mesh;
		using Elem             = typename FunctionSpace::Elem;
		using Dev              = typename FunctionSpace::Device;
		using Vector           = typename FunctionSpace::Vector;
		using Matrix           = typename FunctionSpace::Matrix;
		using Comm             = typename FunctionSpace::Comm;

		static const int Dim    	= Elem::Dim;
		static const int NFunctions = Elem::NFunctions;

		using Point            = typename FunctionSpace::Point;
		using Scalar           = typename FunctionSpace::Scalar;
		using SizeType         = typename FunctionSpace::SizeType;
		using Quadrature       = utopia::Quadrature<Elem, 2>;
		using ElementMatrix    = utopia::StaticMatrix<Scalar, NFunctions, NFunctions>;
		using ElementVector    = utopia::StaticVector<Scalar, NFunctions>;

		using Coefficient = utopia::Coefficient<FunctionSpace>;
		using LEKernel    = utopia::LinearElasticityKernel<Scalar>;

		LinearElasticityFE(FunctionSpace &space, const Scalar &mu = 1.0, const Scalar &lambda = 1.0)
		: space_(space),
		  mu_(mu),
		  lambda_(lambda),
		  debug_(false),
		  x_coeff_(utopia::make_unique<Coefficient>(space_)),
		  grad_(space,  quadrature_),
		  differential_(space, quadrature_)
		{}

		void read(Input &in) override {
			in.get("debug", debug_);
			in.get("mu", mu_);
			in.get("lambda", lambda_);
		}

		inline Size size() const
		{
			const SizeType n_dofs = space_.n_dofs();
			return {n_dofs, n_dofs};
		}

		inline Size local_size() const
		{
			const SizeType n_dofs = space_.n_local_dofs();
			return {n_dofs, n_dofs};
		}

		inline Comm &comm() override { return space_.comm(); }
		inline const Comm &comm() const override { return space_.comm(); }

		bool apply(const Vector &x, Vector &y) const override {			
			const Comm &comm = space_.comm();

			if(y.empty()) {
				space_.create_vector(y);
			} else {
				y.set(0.0);
			}

			Quadrature quadrature;
			auto &&space_view = space_.view_device();

			x_coeff_->update(x);

			LinearElasticity<FunctionSpace, Quadrature> elast(space_, quadrature, mu_, lambda_);

			{
				auto x_view    = x_coeff_->view_device();	
				auto y_view    = space_.assembly_view_device(y);
				auto grad_view = grad_.view_device();
				auto dx_view   = differential_.view_device();

				Dev::parallel_for(
					space_.local_element_range(),
					UTOPIA_LAMBDA(const SizeType &i)
					{
						Elem e;
						ElementVector coeff, el_vec;
						el_vec.set(0.0);

						space_view.elem(i, e);
						x_view.get(e, coeff);

						auto grad = grad_view.make(e);
						auto dx   = dx_view.make(e);
						
						const auto n = grad.n_points();
						for(SizeType k = 0; k < n; ++k) {
						    for(SizeType j = 0; j < grad.n_functions(); ++j) {
						        const auto g_test = grad(j, k);
						        el_vec(j) += LEKernel::apply(
						        	mu_,
						        	lambda_,
						        	coeff(j),
						        	g_test,
						        	g_test,
						        	dx(k)
						        );

						        for(SizeType l = j + 1; l < grad.n_functions(); ++l) {
						        	const auto g_trial = grad(l, k);
						            
						            el_vec(j) += LEKernel::apply(
						            	mu_,
						            	lambda_,
						            	coeff(l),
						            	g_trial,
						            	g_test,
						            	dx(k)
						            );

						            el_vec(l) += LEKernel::apply(
						            	mu_,
						            	lambda_,
						            	coeff(j),
						            	g_test,
						            	g_trial,
						            	dx(k)
						            );
						        }
						    }
						}

						space_view.add_vector(e, el_vec, y_view);
					}
				);
			}

			space_.copy_at_constrained_dofs(x, y);
			return true;
		}

	private:
		FunctionSpace space_;
		Scalar mu_, lambda_;
		bool debug_;

		std::unique_ptr<Coefficient> x_coeff_;

		//assembly
		Quadrature quadrature_;
		PhysicalGradient<FunctionSpace, Quadrature> grad_;
		Differential<FunctionSpace, Quadrature> differential_;
	};

}

#endif //UTOPIA_LINEAR_ELASTICITY_FE_HPP

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
#include "utopia_make_unique.hpp"
#include "utopia_StrainView.hpp"


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
		using Quadrature       = utopia::Quadrature<Elem, 2*(Elem::Order -1)>;
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
		{
			H_ptr_ = utopia::make_unique<Matrix>();
		}

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

			// LinearElasticity<FunctionSpace, Quadrature> elast(space_, quadrature, mu_, lambda_);
			Strain<FunctionSpace, Quadrature> strain(space_, quadrature);

			{
				auto x_view    = x_coeff_->view_device();
				auto y_view    = space_.assembly_view_device(y);
				auto strain_view = strain.view_device();
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

						auto &&strain = strain_view.make(e);
						auto dx       = dx_view.make(e);

						const SizeType n_qp  = strain.n_points();
						const SizeType n_fun = strain.n_functions();

						//taking advantage of symmetry
						for(SizeType qp = 0; qp < n_qp; ++qp) {
							const auto dx_qp = dx(qp);

						    for(SizeType j = 0; j < n_fun; ++j) {
						    	auto &&strain_test  = strain(j, qp);

						    	el_vec(j) += LEKernel::strain_apply(
						    		mu_,
						    		lambda_,
						    		strain_test,
						    		strain_test,
						    		dx_qp
						    	) * coeff(j);

						        for(SizeType l = j + 1; l < n_fun; ++l) {
						        	auto &&strain_trial = strain(l, qp);

						        	const auto v = LEKernel::strain_apply(
						            	mu_,
						            	lambda_,
						            	strain_trial,
						            	strain_test,
						            	dx_qp
						            );

						            el_vec(l) += v * coeff(j);
						            el_vec(j) += v * coeff(l);
						        }
						    }
						}

						space_view.add_vector(e, el_vec, y_view);
					}
				);
			}

			space_.copy_at_constrained_dofs(x, y);
			assert(check_with_matrix(x, y));
			return true;
		}

		bool check_with_matrix(const Vector &x, const Vector &y) const
		{
			if(H_ptr_->empty()) {
				space_.create_matrix(*H_ptr_);
				hessian(x, *H_ptr_);
			}

			Vector mat_y = (*H_ptr_) * x;

			Scalar norm_diff = norm_infty(y - mat_y);

			// std::cout << "norm_diff: " << norm_diff << std::endl;
			assert(norm_diff < 1e-10);
			return norm_diff < 1e-10;
		}

		bool hessian(const Vector &, Matrix &H) const
		{
			LinearElasticity<FunctionSpace, Quadrature> elast(space_, quadrature_, mu_, lambda_);
			MassMatrix<FunctionSpace, Quadrature> mass_matrix(space_, quadrature_);

			{
				auto &&space_view = space_.view_device();
			    auto mat_view     = space_.assembly_view_device(H);
			    auto elast_view   = elast.view_device();

			    if(debug_) {
			        disp("elast");
			        elast_view.describe();
			    }

			    Dev::parallel_for(
			        space_.local_element_range(),
			        UTOPIA_LAMBDA(const SizeType &i)
			    {
			        Elem e;

			        //FIXME this is too big for GPU stack memory for hexas
			        ElementMatrix el_mat;
			        space_view.elem(i, e);

			        //Assemble local elast
			        el_mat.set(0.0);
			        elast_view.assemble(i, e, el_mat);
			        space_view.add_matrix(e, el_mat, mat_view);
			    });
			}

			space_.apply_constraints(H);
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

		std::unique_ptr<Matrix> H_ptr_;
	};

}

#endif //UTOPIA_LINEAR_ELASTICITY_FE_HPP

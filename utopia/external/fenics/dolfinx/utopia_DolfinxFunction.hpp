#ifndef UTOPIA_DOLFINX_FUNCTION_HPP
#define UTOPIA_DOLFINX_FUNCTION_HPP

#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Traits.hpp"

#include "utopia_Function.hpp"

#include <memory>

#include <dolfinx.h>

namespace utopia {

	class DolfinxFunction : public Function<PetscMatrix, PetscVector> {
	public:
		using Scalar = Traits<PetscVector>::Scalar;
		using Form_t = dolfinx::fem::Form<Scalar>;
		using DirichletBC_t = dolfinx::fem::DirichletBC<Scalar>;

		~DolfinxFunction();

		DolfinxFunction(
			const std::shared_ptr<Form_t> &value,
			const std::shared_ptr<Form_t> &gradient,
			const std::shared_ptr<Form_t> &hessian,
			std::vector<std::shared_ptr<const DirichletBC_t>> bcs);

		DolfinxFunction(
			const std::shared_ptr<Form_t> &residual,
			const std::shared_ptr<Form_t> &jacobian,
			std::vector<std::shared_ptr<const DirichletBC_t>> bcs);

		bool hessian(const utopia::PetscVector &x, utopia::PetscMatrix &H) const override;
		bool value(const utopia::PetscVector &x, PetscScalar &value) const override;
		bool gradient(const utopia::PetscVector &x, utopia::PetscVector &g) const override;

	private:
		void init();
		void destory();

		class Impl;
		std::unique_ptr<Impl> impl_;
	};
}

#endif //UTOPIA_DOLFINX_FUNCTION_HPP

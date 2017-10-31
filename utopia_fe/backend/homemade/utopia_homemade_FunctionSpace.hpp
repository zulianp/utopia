#ifndef UTOPIA_HOMEMADE_FUNCTION_SPACE_HPP
#define UTOPIA_HOMEMADE_FUNCTION_SPACE_HPP 

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_homemade_FEForwardDeclarations.hpp"

#include "utopia_FunctionalTraits.hpp"
#include "utopia_FunctionSpace.hpp"
#include "utopia_homemade_Mesh.hpp"
#include "utopia_homemade_AssemblyContext.hpp"

namespace utopia {

	class HMFESpace : public FunctionSpace<HMFESpace> {
	public:
		HMFESpace()
		: mesh_(std::make_shared<Mesh>()) 
		{}

		inline Mesh &mesh()
		{
			return *mesh_;
		}

		inline const Mesh &mesh() const
		{
			return *mesh_;
		}


		std::shared_ptr<Mesh> mesh_;
	};

	template<>
	class Traits<HMFESpace> {
	public:
		static const int Backend = HOMEMADE;
		static const int Order = 1;
		static const int FILL_TYPE = FillType::DENSE;

		typedef double Scalar;
		typedef utopia::Vectord Vector;
		typedef utopia::Matrixd Matrix;

		typedef utopia::HMFESpace Implementation;
		typedef utopia::HMDerivative GradientType;
		typedef utopia::HMFun DivergenceType;
		typedef utopia::HMJacobian JacobianType;
	};

	template<>
	class FunctionalTraits<HMFESpace, AssemblyContext<HOMEMADE> > {
	public:
		inline static int type(const HMFESpace &space,  const AssemblyContext<HOMEMADE> &ctx)
		{
			return utopia::POLYNOMIAL_FUNCTION;
		}

		inline static int order(const HMFESpace &space, const AssemblyContext<HOMEMADE> &ctx)
		{
			return space.mesh().element_order(ctx.current_element);
		}
	};
}

#endif //UTOPIA_HOMEMADE_FUNCTION_SPACE_HPP

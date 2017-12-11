#ifndef UTOPIA_ELASTIC_MATERIAL_HPP
#define UTOPIA_ELASTIC_MATERIAL_HPP 

#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_LameeParameters.hpp"

namespace utopia {

	class ElasticMaterial {
	public:
		typedef utopia::ProductFunctionSpace<LibMeshFunctionSpace> FunctionSpaceT;

		virtual bool init(
			const FunctionSpaceT &V,
			const LameeParameters &params,
			const DVectord &displacement0,
			DSMatrixd &stiffness_matrix,
			DVectord  &internal_stress) = 0;

		virtual bool update(
			const DVectord &displacement,
			DSMatrixd &stiffness_matrix,
			DVectord  &internal_stress) = 0;
	};

	class LinearElasticity : public ElasticMaterial {
	public:
		// typedef ElasticMaterial::FunctionSpaceT FunctionSpaceT;

		virtual bool init(
			const FunctionSpaceT &V,
			const LameeParameters &params,
			const DVectord &displacement0,
			DSMatrixd &stiffness_matrix,
			DVectord  &internal_stress) override;

		bool update(
			const DVectord &displacement,
			DSMatrixd &stiffness_matrix,
			DVectord  &internal_stress) override;


		bool assemble_hessian(
			const FunctionSpaceT &V,
			const LameeParameters &params,
			const DVectord &displacement,
			DSMatrixd &stiffness_matrix);
	};

}

#endif //UTOPIA_ELASTIC_MATERIAL_HPP

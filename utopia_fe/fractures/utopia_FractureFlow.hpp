#ifndef UTOPIA_FRACTURE_FLOW_HPP
#define UTOPIA_FRACTURE_FLOW_HPP

#include "utopia_Input.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"

#include "libmesh/parallel_mesh.h"

namespace utopia {

	class FractureFlow final : public Configurable {
	public:
		typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
		typedef utopia::Traits<FunctionSpaceT> TraitsT;
		typedef typename TraitsT::Matrix ElementMatrix;


		FractureFlow(libMesh::Parallel::Communicator &comm);
		void read(Input &is) override;
		bool empty() const;
		void describe(std::ostream &os = std::cout) const;

		UIMesh<libMesh::DistributedMesh> mesh;
		UIFunctionSpace<FunctionSpaceT>  space;
		std::shared_ptr< UIForcingFunction<FunctionSpaceT, UVector> > forcing_function;
		std::shared_ptr<UIFunction<double>> sampler;

		ElementMatrix diffusion_tensor;
	};

	class FractureFlowAuxSystem {
	public:
	    FractureFlowAuxSystem(LibMeshFunctionSpace &V, const std::string &name = "aperture");
	    void sample(const std::shared_ptr<UIFunction<double>> &sampler);

	    UVector sampled;

	private:

	    libMesh::LinearImplicitSystem &aux_;
	    std::vector<int> var_nums_;
	};

}

#endif //UTOPIA_FRACTURE_FLOW_HPP

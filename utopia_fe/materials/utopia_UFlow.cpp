#include "utopia_UFlow.hpp"

#include "utopia_Flow.hpp"
#include "utopia_NewFlow.hpp"
#include "utopia_FlowWithFractures.hpp"

namespace utopia {

	template<class FunctionSpace, class Matrix, class Vector>
	UFlow<FunctionSpace, Matrix, Vector>::UFlow(FunctionSpace &space)
	: space_(space), rescale_(1.0)
	{}

	template<class FunctionSpace, class Matrix, class Vector>
	void UFlow<FunctionSpace, Matrix, Vector>::read(Input &in)
	{
		//FIXME allow also other models
		std::string type;
		in.get("type", type);
		if("FlowWithFractures" == type)  {
		    auto flow = utopia::make_unique<FlowWithFractures<FunctionSpace, Matrix, Vector> >(space_);
		    flow->rescale(rescale_);
		    flow->read(in);
		    flow_model_ = std::move(flow);
		} else if("NewFlow" == type) {
		    auto flow = utopia::make_unique<NewFlow<FunctionSpace, Matrix, Vector> >(space_);
		    flow->rescale(rescale_);
		    flow->read(in);
		    flow_model_ = std::move(flow);
		} else {
		    auto flow = utopia::make_unique<Flow<FunctionSpace, Matrix, Vector> >(space_);
		    flow->rescale(rescale_);
		    flow->read(in);
		    flow_model_ = std::move(flow);
		}
	}

	template<class FunctionSpace, class Matrix, class Vector>
	bool UFlow<FunctionSpace, Matrix, Vector>::assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient)
	{
		assert(is_valid());
		return flow_model_->assemble_hessian_and_gradient(x, hessian, gradient);
	}

	template class UFlow<LibMeshFunctionSpace, USparseMatrix, UVector>;
}

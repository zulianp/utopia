#ifndef UTOPIA_EMBEDDED_MODEL_HPP
#define UTOPIA_EMBEDDED_MODEL_HPP

#include "utopia_Model.hpp"
#include "utopia_FEModel.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_Flow.hpp"

#include "libmesh/parallel_mesh.h"

#include <memory>

namespace utopia {

    template<class Matrix, class Vector>
    class DiscreteFractureNetwork final : /*public FEModel<LibMeshFunctionSpace, Matrix, Vector>*/ public Configurable {
    public:
        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        void read(Input &in) override
        {
            in.get("mesh", mesh_);
            in.get("space", space_);

            //FIXME
            flow_model_ = std::make_shared<Flow<FunctionSpaceT, Matrix, Vector> >(space_.space().subspace(0));
            flow_model_->read(in);
        }

        inline bool assemble_flow(const Vector &x, Matrix &hessian, Vector &gradient)
        {
            assert(flow_model_);
            
            bool ok = flow_model_->assemble_hessian_and_gradient(x, hessian, gradient);

            //TODO handle non-conformities

            return ok;
        }

        DiscreteFractureNetwork(libMesh::Parallel::Communicator &comm)
        : mesh_(comm), space_(make_ref(mesh_))
        {}

        inline FunctionSpaceT &space()
        {
            return space_.space().subspace(0);
        }

        inline const FunctionSpaceT &space() const
        {
            return space_.space().subspace(0);
        }

    private:

        UIMesh<libMesh::DistributedMesh> mesh_;
        UIFunctionSpace<FunctionSpaceT>  space_;
        std::shared_ptr<Model<Matrix, Vector>> flow_model_;
        //TODO add Mortar here
    };

}

#endif //UTOPIA_EMBEDDED_MODEL_HPP

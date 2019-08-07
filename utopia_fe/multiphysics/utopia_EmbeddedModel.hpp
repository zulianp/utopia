#ifndef UTOPIA_EMBEDDED_MODEL_HPP
#define UTOPIA_EMBEDDED_MODEL_HPP

#include "utopia_Model.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_Flow.hpp"

#include "libmesh/parallel_mesh.h"

#include <memory>

namespace utopia {

    template<class Matrix, class Vector>
    class EmbeddedModel final : public Model<Matrix, Vector> {
    public:
        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        void read(Input &in) override
        {
            in.get("mesh", mesh_);
            in.get("space", space_);

            //FIXME
            model_ = std::make_shared<Flow<FunctionSpaceT, Matrix, Vector> >(space_.space().subspace(0));

        }

        inline bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        {
            assert(model_);
            
            bool ok = model_->assemble_hessian_and_gradient(x, hessian, gradient);

            //TODO handle non-conformities

            return ok;
        }

        UIMesh<libMesh::DistributedMesh> mesh_;
        UIFunctionSpace<FunctionSpaceT>  space_;
        std::shared_ptr<Model<Matrix, Vector>> model_;
        //TODO add Mortar here
    };

}

#endif //UTOPIA_EMBEDDED_MODEL_HPP

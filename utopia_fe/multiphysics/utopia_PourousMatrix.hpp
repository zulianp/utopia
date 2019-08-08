#ifndef UTOPIA_BACKGROUND_MODEL_HPP
#define UTOPIA_BACKGROUND_MODEL_HPP

#include "utopia_FEModel.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_Flow.hpp"

#include "libmesh/parallel_mesh.h"

#include <memory>

namespace utopia {

    template<class Matrix, class Vector>
    class PourousMatrix final : /*public FEModel<LibMeshFunctionSpace, Matrix, Vector> */ public Configurable  {
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

        inline bool init()
        {
            // NewTransferAssembler new_assembler;
            // new_assembler.use_convert_transfer(params_->use_convert_transfer);
            // new_assembler.remove_incomplete_intersections(params_->remove_incomplete_intersections);

            // if(params_->surface_transfer) {
            //     if(!new_assembler.surface_assemble(from_mesh, from_dofs, to_mesh, to_dofs, opts)) {
            //         return false;
            //     }
            // } else {
            //     if(!new_assembler.assemble(from_mesh, from_dofs, to_mesh, to_dofs, opts)) {
            //         return false;
            //     }
            // }

            // operator_ = new_assembler.build_operator();
        }

        inline FunctionSpaceT &space()
        {
            return space_.space().subspace(0);
        }

        inline const FunctionSpaceT &space() const
        {
            return space_.space().subspace(0);
        }

        PourousMatrix(libMesh::Parallel::Communicator &comm)
        : mesh_(comm), space_(make_ref(mesh_))
        {}

    private:
        UIMesh<libMesh::DistributedMesh> mesh_;
        UIFunctionSpace<FunctionSpaceT>  space_;
        std::shared_ptr<Model<Matrix, Vector>> flow_model_;
        //TODO add Mortar here
    };

}

#endif //UTOPIA_BACKGROUND_MODEL_HPP

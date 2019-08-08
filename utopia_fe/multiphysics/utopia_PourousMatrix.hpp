#ifndef UTOPIA_BACKGROUND_MODEL_HPP
#define UTOPIA_BACKGROUND_MODEL_HPP

#include "utopia_FEModel.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_Flow.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_NewTransferAssembler.hpp"

#include "libmesh/parallel_mesh.h"

#include <memory>

namespace utopia {

    template<class Matrix, class Vector>
    class Mortar : public Configurable {
    public:

        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        using Scalar   = UTOPIA_SCALAR(Vector);
        using SizeType = UTOPIA_SIZE_TYPE(Vector);

        void read(Input &in) override
        {
            opts.read(in);
        }

        bool init(FunctionSpaceT &space)
        {
            empty_ = false;
            if(opts.tags.empty()) {
                empty_ = true;
                return false;
            }

            NewTransferAssembler mortar_assembler;
            if(!mortar_assembler.surface_assemble(space.mesh(), space.dof_map(), opts)) {
                empty_ = true;
                std::cerr << "[Error] failed to find any surface to couple" << std::endl;
                return false;
            }

            auto op = mortar_assembler.build_operator();
            mortar_matrix_ = op->matrix();

            is_constrained_ = sum(*mortar_matrix_, 1);
            is_unconstrained_ = is_constrained_;

            each_transform(
                is_unconstrained_,
                is_unconstrained_, 
                [](const SizeType i, const Scalar val) -> Scalar {
                    if(val > 0.9) {
                        return 0.0;
                    } else {
                        return 1.0;
                    }
                }
            );

            // disp(is_constrained_);
            // disp(is_unconstrained_);
            return true;
        }

        void constrain_system(Matrix &A, Vector &b)
        {
            const auto &T = *mortar_matrix();

            Matrix slave_A = transpose(T) * A * T;
            Vector slave_b = transpose(T) * b;
            
            A += slave_A;
            b += slave_b;

            set_zero_rows(A, is_constrained_, 1.0);
            b = e_mul(is_unconstrained_, b);
        }

        void unconstrain_solution(Vector &x) {
            const auto &T = *mortar_matrix();
            x += T * x;

            x = is_constrained_;
        }

        inline std::shared_ptr<Matrix> mortar_matrix()
        {
            return mortar_matrix_;
        }

        inline std::shared_ptr<const Matrix> mortar_matrix() const
        {
            return mortar_matrix_;
        }

        inline bool empty() const
        {
            return empty_;
        }

    private:
        TransferOptions opts;
        std::shared_ptr<Matrix> mortar_matrix_;
        Vector is_constrained_, is_unconstrained_;
        bool empty_;
    };


    template<class Matrix, class Vector>
    class PourousMatrix final : /*public FEModel<LibMeshFunctionSpace, Matrix, Vector> */ public Configurable  {
    public:
        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        using Mortar = utopia::Mortar<Matrix, Vector>;

        void read(Input &in) override
        {
            in.get("mesh", mesh_);
            in.get("space", space_);

            //FIXME allow also other models
            flow_model_ = std::make_shared<Flow<FunctionSpaceT, Matrix, Vector> >(space_.space().subspace(0));
            flow_model_->read(in);

            in.get("mortar", mortar_);

            init();
        }

        inline bool assemble_flow(const Vector &x, Matrix &hessian, Vector &gradient)
        {
            assert(flow_model_);

            if(!flow_model_->assemble_hessian_and_gradient(x, hessian, gradient)) return false;
            if(mortar_.empty()) return true;
            mortar_.constrain_system(hessian, gradient);
            return true;
        }

        inline void disassemble_flow(Vector &x)
        {
            if(mortar_.empty()) return;
            mortar_.unconstrain_solution(x);
        }

        inline void init()
        {
            mortar_.init(space());
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
        Mortar mortar_;
        //TODO add Mortar here
    };

}

#endif //UTOPIA_BACKGROUND_MODEL_HPP

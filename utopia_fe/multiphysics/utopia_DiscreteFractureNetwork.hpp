#ifndef UTOPIA_EMBEDDED_MODEL_HPP
#define UTOPIA_EMBEDDED_MODEL_HPP

#include "utopia_Model.hpp"
#include "utopia_FEModel.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_Flow.hpp"
#include "utopia_LowerDimTransfer.hpp"

#include "libmesh/parallel_mesh.h"

#include <memory>

namespace utopia {


    template<class Matrix, class Vector>
    class FractureCoupling final : public Configurable {
    public:

        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        using Scalar   = UTOPIA_SCALAR(Vector);
        using SizeType = UTOPIA_SIZE_TYPE(Vector);

        void read(Input &in) override
        {
            opts.read(in);
            in.get("export-constrained", export_constrained_);
        }

        bool init(FunctionSpaceT &space)
        {
            empty_ = false;
            
            if(opts.tags.empty()) {
                empty_ = true;
                return false;
            }

            LowerDimTransfer transfer_assembler;
            if(!transfer_assembler.assemble(space.mesh(), space.dof_map(), opts)) {
                empty_ = true;
                std::cerr << "[Error] failed to find any surface to couple" << std::endl;
                return false;
            }

            auto op = transfer_assembler.build_operator();
            coupling_matrix_ = op->matrix();

            is_constrained_ = sum(*coupling_matrix_, 1);
            is_unconstrained_ = is_constrained_;

            //create a clean boolean vector

            each_transform(
                is_constrained_,
                is_constrained_, 
                [](const SizeType i, const Scalar val) -> Scalar {
                    if(val > 0.99) {
                        return 1.0;
                    } else {
                        return 0.0;
                    }
                }
            );

            each_transform(
                is_unconstrained_,
                is_unconstrained_, 
                [](const SizeType i, const Scalar val) -> Scalar {
                    if(val > 0.99) {
                        return 0.0;
                    } else {
                        return 1.0;
                    }
                }
            );

            if(export_constrained_) {
                write("constrained.e",   space, is_constrained_);
                write("unconstrained.e", space, is_unconstrained_);
                write("T.m", *coupling_matrix_);
            }

            *coupling_matrix_ += local_identity(local_size(*coupling_matrix_));

            return true;
        }

        void constrain_system(Matrix &A, Vector &b)
        {
            const auto &T = *coupling_matrix();

            A = transpose(T) * A * T;
            b = transpose(T) * b;
            
            set_zero_rows(A, is_constrained_, 1.0);
            b = e_mul(is_unconstrained_, b);
        }

        void unconstrain_solution(Vector &x) {
            const auto &T = *coupling_matrix();
            x = T * x;
        }

        inline std::shared_ptr<Matrix> coupling_matrix()
        {
            return coupling_matrix_;
        }

        inline std::shared_ptr<const Matrix> coupling_matrix() const
        {
            return coupling_matrix_;
        }

        inline bool empty() const
        {
            return empty_;
        }

        FractureCoupling()
        : empty_(true), export_constrained_(false)
        {}

    private:
        TransferOptions opts;
        std::shared_ptr<Matrix> coupling_matrix_;
        Vector is_constrained_, is_unconstrained_;
        bool empty_;
        bool export_constrained_;
    };


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

            in.get("coupling", coupling_);

            init();
        }

        inline bool assemble_flow(const Vector &x, Matrix &hessian, Vector &gradient)
        {
            assert(flow_model_);
            
            if(!flow_model_->assemble_hessian_and_gradient(x, hessian, gradient)) return false;

            if(coupling_.empty()) return true;

            coupling_.constrain_system(hessian, gradient);
            return true;
        }

        inline void disassemble_flow(Vector &x)
        {
            if(coupling_.empty()) return;

            coupling_.unconstrain_solution(x);
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

        inline void init()
        {
            coupling_.init(space());
        }
        
    private:

        UIMesh<libMesh::DistributedMesh> mesh_;
        UIFunctionSpace<FunctionSpaceT>  space_;
        std::shared_ptr<Model<Matrix, Vector>> flow_model_;
        FractureCoupling<Matrix, Vector> coupling_;
    };

}

#endif //UTOPIA_EMBEDDED_MODEL_HPP

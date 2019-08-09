#ifndef UTOPIA_BACKGROUND_MODEL_HPP
#define UTOPIA_BACKGROUND_MODEL_HPP

#include "utopia_FEModel.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_Flow.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_NewTransferAssembler.hpp"
#include "utopia_FluxPostProcessor.hpp"

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
            in.get("export-constrained", export_constrained_);
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

            *mortar_matrix_ += local_identity(local_size(*mortar_matrix_));

            if(export_constrained_) {
                write("constrained.e",   space, is_constrained_);
                write("unconstrained.e", space, is_unconstrained_);
                write("T.m", *mortar_matrix_);
            }

            return true;
        }

        void constrain_system(Matrix &A, Vector &b)
        {
            const auto &T = *mortar_matrix();

            A = transpose(T) * A * T;
            b = transpose(T) * b;
            
            set_zero_rows(A, is_constrained_, 1.0);
            b = e_mul(is_unconstrained_, b);
        }

        void unconstrain_solution(Vector &x) {
            const auto &T = *mortar_matrix();
            x = T * x;
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

        Mortar()
        : export_constrained_(false)
        {}

    private:
        TransferOptions opts;
        std::shared_ptr<Matrix> mortar_matrix_;
        Vector is_constrained_, is_unconstrained_;
        bool empty_;
        bool export_constrained_;
    };

    template<class Matrix, class Vector>
    class PostProcessable : public virtual Configurable {
    public:
        using FunctionSpaceT = utopia::LibMeshFunctionSpace;
        using PostProcessor  = utopia::PostProcessor<FunctionSpaceT, Vector>;

        virtual ~PostProcessable() {}

        virtual void read(Input &in) override {

            in.get("post-processors", [this](Input &in) {
                in.get_all([this](Input &in) {
                    std::string type;
                    in.get("type", type);

                    if(type == "flux") {
                        auto flux = std::make_shared<FluxPostProcessor<FunctionSpaceT, UVector>>();    
                        flux->read(in);

                        post_processors_.push_back(flux);

                    } else if(type == "avg") {
                        auto flux = std::make_shared<AverageHeadPostProcessor<FunctionSpaceT, UVector>>();              
                        flux->read(in);

                        post_processors_.push_back(flux);
                    }
                });
            });

        }

        virtual void post_process(FunctionSpaceT &space, const Vector &x)
        {
            for(auto pp : post_processors_) {
                pp->apply(space, x);
            }
        }

    private:
        std::vector<std::shared_ptr<PostProcessor>> post_processors_;
    };


    template<class Matrix, class Vector>
    class PourousMatrix final : /*public FEModel<LibMeshFunctionSpace, Matrix, Vector> */ public PostProcessable<Matrix, Vector>  {
    public:
        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        using Mortar = utopia::Mortar<Matrix, Vector>;
        using Super  = utopia::PostProcessable<Matrix, Vector>;
        

        void read(Input &in) override
        {
            Super::read(in);

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

        inline void post_process_flow(const Vector &x)
        {
            Super::post_process(space(), x);
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
        
    };

}

#endif //UTOPIA_BACKGROUND_MODEL_HPP

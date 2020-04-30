#ifndef UTOPIA_BACKGROUND_MODEL_HPP
#define UTOPIA_BACKGROUND_MODEL_HPP

#include "utopia_FEModel.hpp"
#include "utopia_FluxPostProcessor.hpp"
#include "utopia_FractureFlowUtils.hpp"
#include "utopia_LinePostProcessor.hpp"
#include "utopia_NewTransferAssembler.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_UFlow.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"

#include "libmesh/parallel_mesh.h"

#include <memory>

namespace utopia {

    template <class Matrix, class Vector>
    class Mortar : public Configurable {
    public:
        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        using Scalar = UTOPIA_SCALAR(Vector);
        using SizeType = UTOPIA_SIZE_TYPE(Vector);

        void read(Input &in) override {
            opts.read(in);
            in.get("export-constrained", export_constrained_);
        }

        bool init(FunctionSpaceT &space) {
            empty_ = false;
            if (opts.tags.empty()) {
                empty_ = true;
                return false;
            }

            NewTransferAssembler mortar_assembler;
            if (!mortar_assembler.surface_assemble(space.mesh(), space.dof_map(), opts)) {
                empty_ = true;
                std::cerr << "[Error] failed to find any surface to couple" << std::endl;
                return false;
            }

            auto op = mortar_assembler.build_operator();
            transfer_matrix_ = op->matrix();

            is_constrained_ = sum(*transfer_matrix_, 1);
            is_unconstrained_ = is_constrained_;

            // create a clean boolean vector

            each_transform(is_constrained_, is_constrained_, [](const SizeType i, const Scalar val) -> Scalar {
                if (val > 0.99) {
                    return 1.0;
                } else {
                    return 0.0;
                }
            });

            each_transform(is_unconstrained_, is_unconstrained_, [](const SizeType i, const Scalar val) -> Scalar {
                if (val > 0.99) {
                    return 0.0;
                } else {
                    return 1.0;
                }
            });

            if (export_constrained_) {
                write("constrained.e", space, is_constrained_);
                write("unconstrained.e", space, is_unconstrained_);

                rename("m", *transfer_matrix_);
                write("M.m", *transfer_matrix_);
            }

            set_zero_at_constraint_rows(space.dof_map(), *transfer_matrix_);

            mortar_matrix_ = std::make_shared<USparseMatrix>();
            *mortar_matrix_ = *transfer_matrix_ + local_identity(local_size(*transfer_matrix_));
            return true;
        }

        void constrain_system(Matrix &A, Vector &b) {
            assert(mortar_matrix());
            const auto &T = *mortar_matrix();

            A = transpose(T) * A * T;
            b = transpose(T) * b;

            set_zero_rows(A, is_constrained_, 1.0);
            b = e_mul(is_unconstrained_, b);
        }

        void unconstrain_solution(Vector &x) {
            assert(mortar_matrix());
            const auto &T = *mortar_matrix();
            x = T * x;
        }

        inline std::shared_ptr<Matrix> mortar_matrix() { return mortar_matrix_; }

        inline std::shared_ptr<const Matrix> mortar_matrix() const { return mortar_matrix_; }

        inline std::shared_ptr<const Matrix> transfer_matrix() const { return transfer_matrix_; }

        inline bool empty() const { return empty_; }

        inline const Vector &is_constrained() const { return is_constrained_; }

        inline void compute_mortar_matrix_without_slave_dofs() {
            mortar_matrix_without_slave_dofs_ = std::make_shared<USparseMatrix>();
            *mortar_matrix_without_slave_dofs_ = *transfer_matrix_ + USparseMatrix(diag(is_unconstrained_));
        }

        inline std::shared_ptr<Matrix> mortar_matrix_without_slave_dofs() { return mortar_matrix_without_slave_dofs_; }

        inline std::shared_ptr<const Matrix> mortar_matrix_without_slave_dofs() const {
            return mortar_matrix_without_slave_dofs_;
        }

        Mortar() : empty_(true), export_constrained_(false) {}

    private:
        TransferOptions opts;
        std::shared_ptr<Matrix> mortar_matrix_;
        std::shared_ptr<Matrix> transfer_matrix_;
        std::shared_ptr<Matrix> mortar_matrix_without_slave_dofs_;
        Vector is_constrained_, is_unconstrained_;
        bool empty_;
        bool export_constrained_;
    };

    template <class Matrix, class Vector>
    class GradientRecovery : public Configurable {
    public:
        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        using Scalar = UTOPIA_SCALAR(Vector);

        GradientRecovery();

        void read(Input &in) override;

        inline bool empty() const { return grad_space_.empty(); }

        inline int n_refinements() const { return n_refinements_; }
        /// calls estimate_error then apply_refinement
        bool refine(const Mortar<Matrix, Vector> &mortar, FunctionSpaceT &V, const Vector &sol);
        void estimate_error(const Mortar<Matrix, Vector> &mortar, FunctionSpaceT &V, const Vector &sol);
        bool apply_refinement(FunctionSpaceT &V);

    private:
        UVector all_values_;
        UVector error_;
        ProductFunctionSpace<FunctionSpaceT> grad_space_;
        int system_num_;
        int error_var_num_;
        Scalar max_local_error_;
        int n_refinements_;
        int max_refinements_;

        void init(FunctionSpaceT &V);
        void append_error_estimate(FunctionSpaceT &V);
    };

    template <class Matrix, class Vector>
    class PostProcessable : public virtual Configurable {
    public:
        using FunctionSpaceT = utopia::LibMeshFunctionSpace;
        using PostProcessor = utopia::PostProcessor<FunctionSpaceT, Vector>;
        using Scalar = UTOPIA_SCALAR(Vector);

        virtual ~PostProcessable() {}

        virtual void read(Input &in) override {
            in.get("post-processors", [this](Input &in) {
                in.get_all([this](Input &in) {
                    std::string type;
                    in.get("type", type);

                    std::cout << "post-processor type: " << type << std::endl;

                    if (type == "flux") {
                        auto flux = std::make_shared<FluxPostProcessor<FunctionSpaceT, UVector>>();
                        // flux->rescale(rescale());
                        flux->read(in);

                        post_processors_.push_back(flux);

                    } else if (type == "avg") {
                        auto flux = std::make_shared<AverageHeadPostProcessor<FunctionSpaceT, UVector>>();
                        // flux->rescale(rescale());
                        flux->read(in);

                        post_processors_.push_back(flux);
                    } else if (type == "sample-line") {
                        auto line_pp = std::make_shared<LinePostProcessor<FunctionSpaceT, UVector>>();
                        // line_pp->rescale(rescale());
                        line_pp->read(in);

                        post_processors_.push_back(line_pp);
                    }
                });
            });
        }

        virtual void post_process(FunctionSpaceT &space, const Vector &x) {
            for (auto pp : post_processors_) {
                pp->apply(space, x);
            }
        }

        inline void rescale(const Scalar rescale) { rescale_ = rescale; }

        inline Scalar rescale() const { return rescale_; }

        PostProcessable() : rescale_(1.0) {}

    private:
        std::vector<std::shared_ptr<PostProcessor>> post_processors_;
        Scalar rescale_;
    };

    template <class Matrix, class Vector>
    class PourousMatrix final
        : /*public FEModel<LibMeshFunctionSpace, Matrix, Vector> */ public PostProcessable<Matrix, Vector> {
    public:
        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        using Mortar = utopia::Mortar<Matrix, Vector>;
        using Super = utopia::PostProcessable<Matrix, Vector>;
        using GradientRecoveryT = utopia::GradientRecovery<Matrix, Vector>;
        using Scalar = UTOPIA_STORE(Vector);

        void read(Input &in) override {
            Super::read(in);

            in.get("mesh", mesh_);
            in.get("space", space_);

            auto flow = std::make_shared<UFlow<FunctionSpaceT, Matrix, Vector>>(space_.space().subspace(0));
            flow->rescale(this->rescale());
            flow->read(in);
            flow_model_ = flow;

            in.get("mortar", mortar_);
            in.get("gradient-recovery", gradient_recovery_);

            init();
        }

        inline bool write(Vector &x) {
            if (gradient_recovery_.empty()) {
                utopia::write(space().equation_system().name() + ".e", space(), x);
            } else {
                utopia::write(
                    space().equation_system().name() + "-" + std::to_string(gradient_recovery_.n_refinements()) + ".e",
                    space(),
                    x);
            }

            return true;
        }

        inline bool refine(const Vector &x) {
            gradient_recovery_.estimate_error(mortar_, space(), x);

            if (gradient_recovery_.apply_refinement(space())) {
                init();
                return true;
            }

            return false;
        }

        inline bool assemble_flow(const Vector &x, Matrix &hessian, Vector &gradient) {
            assert(flow_model_);

            if (!flow_model_->assemble_hessian_and_gradient(x, hessian, gradient)) return false;
            if (mortar_.empty()) return true;

            mortar_.constrain_system(hessian, gradient);
            return true;
        }

        inline void disassemble_flow(Vector &x) {
            if (mortar_.empty()) return;

            mortar_.unconstrain_solution(x);
        }

        inline void init() { mortar_.init(space()); }

        inline void post_process_flow(const Vector &x) { Super::post_process(space(), x); }

        inline FunctionSpaceT &space() { return space_.space().subspace(0); }

        inline const FunctionSpaceT &space() const { return space_.space().subspace(0); }

        inline bool has_mortar_constraints() const { return !mortar_.empty(); }

        Mortar &mortar() { return mortar_; }

        PourousMatrix(libMesh::Parallel::Communicator &comm) : mesh_(comm), space_(make_ref(mesh_)) {}

    private:
        UIMesh<libMesh::DistributedMesh> mesh_;
        UIFunctionSpace<FunctionSpaceT> space_;
        std::shared_ptr<Model<Matrix, Vector>> flow_model_;
        Mortar mortar_;
        GradientRecoveryT gradient_recovery_;
    };

}  // namespace utopia

#endif  // UTOPIA_BACKGROUND_MODEL_HPP

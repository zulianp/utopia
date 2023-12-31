#ifndef UTOPIA_EMBEDDED_MODEL_HPP
#define UTOPIA_EMBEDDED_MODEL_HPP

#include "utopia_FEModel.hpp"
#include "utopia_LowerDimTransfer.hpp"
#include "utopia_Model.hpp"
#include "utopia_UFlow.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"

#include "libmesh/parallel_mesh.h"

#include <memory>

namespace utopia {

    template <class Matrix, class Vector>
    class FractureCoupling final : public Configurable {
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

            LowerDimTransfer transfer_assembler;
            if (!transfer_assembler.assemble(space.mesh(), space.dof_map(), opts)) {
                empty_ = true;
                std::cerr << "[Error] failed to find any surface to couple" << std::endl;
                return false;
            }

            auto op = transfer_assembler.build_operator();
            coupling_matrix_ = op->matrix();

            is_constrained_ = sum(*coupling_matrix_, 1);
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
                write("T.m", *coupling_matrix_);
            }

            *coupling_matrix_ += local_identity(local_size(*coupling_matrix_));

            return true;
        }

        void constrain_system(Matrix &A, Vector &b, const bool identity_on_constrained_dofs = false) {
            const auto &T = *coupling_matrix();

            A = transpose(T) * A * T;
            b = transpose(T) * b;

            if (identity_on_constrained_dofs) {
                set_zero_rows(A, is_constrained_, 1.0);
            } else {
                set_zero_rows(A, is_constrained_, 0.0);
            }

            b = e_mul(is_unconstrained_, b);
        }

        void unconstrain_solution(Vector &x) {
            const auto &T = *coupling_matrix();
            x = T * x;
        }

        inline std::shared_ptr<Matrix> coupling_matrix() { return coupling_matrix_; }

        inline std::shared_ptr<const Matrix> coupling_matrix() const { return coupling_matrix_; }

        inline bool empty() const { return empty_; }

        inline bool refine(const Vector &) { return false; }

        FractureCoupling() : empty_(true), export_constrained_(false) {}

    private:
        TransferOptions opts;
        std::shared_ptr<Matrix> coupling_matrix_;
        Vector is_constrained_, is_unconstrained_;
        bool empty_;
        bool export_constrained_;
    };

    template <class Matrix, class Vector>
    class DiscreteFractureNetwork final : /*public FEModel<LibMeshFunctionSpace, Matrix, Vector>*/ public Configurable {
    public:
        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;
        using Scalar = UTOPIA_SCALAR(Vector);

        void read(Input &in) override {
            in.get("mesh", mesh_);
            in.get("space", space_);
            in.get("identity-on-constrained-dofs", identity_on_constrained_dofs_);

            auto flow = std::make_shared<UFlow<FunctionSpaceT, Matrix, Vector>>(space_.space().subspace(0));
            flow->rescale(rescale_);
            flow->read(in);
            flow_model_ = flow;

            in.get("intersection", coupling_);

            init();
        }

        inline void residual(const Matrix &A, const Vector &b, const Vector &x, Vector &r) const { r = b - A * x; }

        inline bool assemble_flow(const Vector &x, Matrix &hessian, Vector &gradient) {
            assert(flow_model_);

            if (!flow_model_->assemble_hessian_and_gradient(x, hessian, gradient)) return false;

            if (coupling_.empty()) return true;

            coupling_.constrain_system(hessian, gradient, identity_on_constrained_dofs_);
            return true;
        }

        inline void disassemble_flow(Vector &x) {
            if (coupling_.empty()) return;

            coupling_.unconstrain_solution(x);
        }

        DiscreteFractureNetwork(libMesh::Parallel::Communicator &comm)
            : mesh_(comm), space_(make_ref(mesh_)), identity_on_constrained_dofs_(false), rescale_(1) {}

        inline FunctionSpaceT &space() { return space_.space().subspace(0); }

        inline const FunctionSpaceT &space() const { return space_.space().subspace(0); }

        inline void init() { coupling_.init(space()); }

        bool compute_flow() {
            auto &dof_map = space().dof_map();

            Matrix A;
            Vector rhs;
            Vector x = local_zeros(dof_map.n_local_dofs());

            if (!assemble_flow(x, A, rhs)) {
                std::cerr << "[Error] failed to assemble" << std::endl;
                return false;
            }

            apply_boundary_conditions(dof_map, A, rhs);

            Factorization<Matrix, Vector> solver;

            if (!solver.solve(A, rhs, x)) {
                std::cerr << "[Error] failed to solve" << std::endl;
                write("A.m", A);
                write("b.m", rhs);
                return false;
            }

            disassemble_flow(x);
            return export_flow(x);
        }

        inline bool export_flow(Vector &x) {
            const std::string name = space().equation_system().name();
            write(name + ".e", space(), x);
            return true;
        }

        inline void rescale(const Scalar rescale) { rescale_ = rescale; }

    private:
        UIMesh<libMesh::DistributedMesh> mesh_;
        UIFunctionSpace<FunctionSpaceT> space_;
        std::shared_ptr<Model<Matrix, Vector>> flow_model_;
        FractureCoupling<Matrix, Vector> coupling_;
        bool identity_on_constrained_dofs_;
        Scalar rescale_;
    };

}  // namespace utopia

#endif  // UTOPIA_EMBEDDED_MODEL_HPP

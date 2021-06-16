#ifndef UTOPIA_NC_FUNCTION_SPACE_HPP
#define UTOPIA_NC_FUNCTION_SPACE_HPP

#include "utopia_fe_base.hpp"

#include "utopia_IPTransfer.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Transfer.hpp"

#include "utopia_FunctionSpaceBase.hpp"
#include "utopia_fe_Core.hpp"

#ifdef UTOPIA_WITH_MOONOLITH
#include "utopia_FETransferOptions.hpp"
#endif  // UTOPIA_WITH_MOONOLITH

#include <memory>
#include <string>

namespace utopia {

    template <class FunctionSpace>
    class Traits<NCFunctionSpace<FunctionSpace>> : public Traits<FunctionSpace> {};

    template <class FunctionSpace>
    class NCFunctionSpace : public FunctionSpaceBase<typename Traits<FunctionSpace>::Mesh> {
    public:
        using Communicator = typename Traits<FunctionSpace>::Communicator;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using IndexArray = typename Traits<FunctionSpace>::IndexArray;
        using SizeType = typename Traits<FunctionSpace>::SizeType;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using IPTransfer = utopia::IPTransfer<Matrix, Vector>;
        using Mesh = typename Traits<FunctionSpace>::Mesh;

        std::shared_ptr<FunctionSpace> unconstrained_space() const { return space_; }

        NCFunctionSpace(const std::shared_ptr<FunctionSpace> &space) : space_(space) {}
        NCFunctionSpace(const Communicator &comm) : space_(std::make_shared<FunctionSpace>(comm)) {}

        void read(Input &in) override {
            if (!space_) {
                space_ = std::make_shared<FunctionSpace>();
                space_->read(in);
            } else if (space_->empty()) {
                space_->read(in);
            }

            std::string type;
            in.require("type", type);

#ifdef UTOPIA_WITH_MOONOLITH
            if (type == "mortar") {
                FETransfer<FunctionSpace> transfer;
                in.get("mortar", transfer);
                if (!transfer.init_from_decomposition(space_)) {
                    Utopia::Abort("NCFunctionSpace: failed to initialize mortar operator!");
                }

                auto temp = transfer.template build_transfer<IPTransfer>();

                space_->apply_constraints(*temp->I_ptr(), 0.0);

                Vector is_constrained = sum(temp->I(), 1);
                constrained_indices_.reserve(is_constrained.local_size());
                is_constrained.read([this](const SizeType &i, const Scalar &v) {
                    if (v > 0.99) {
                        constrained_indices_.push_back(i);
                    }
                });

                projector_ = temp;
            }
#endif  // UTOPIA_WITH_MOONOLITH
        }

        void restrict(const Matrix &in, Matrix &out) const {
            if (projector_) {
                Matrix temp;
                projector_->restrict(in, temp);
                out += temp;
            }
        }

        void restrict(const Vector &in, Vector &out) const {
            if (projector_) {
                projector_->restrict(in, out);
            }
        }

        // void interpolate(const Matrix &in, Matrix &out) const {
        //     out = in;
        //     if (projector_) {
        //         Vector temp;
        //         projector_->apply(in, temp);
        //         out += temp;
        //     }
        // }

        void interpolate(const Vector &in, Vector &out) const {
            out = in;
            if (projector_) {
                Vector temp;
                projector_->interpolate(in, temp);
                out += temp;
            }
        }

        void project(const Vector &in, Vector &out) const {
            if (projector_) {
                projector_->interpolate(in, out);
            }
        }

        void init(const std::shared_ptr<Mesh> &mesh) override { unconstrained_space()->init(mesh); }
        bool write(const Path &path, const Vector &x) override { return unconstrained_space()->write(path, x); }

        std::shared_ptr<Mesh> mesh_ptr() const override { return unconstrained_space()->mesh_ptr(); }
        const Mesh &mesh() const override { return unconstrained_space()->mesh(); }
        Mesh &mesh() override { return unconstrained_space()->mesh(); }

        const Communicator &comm() const override { return unconstrained_space()->comm(); }

        SizeType n_dofs() const override { return unconstrained_space()->n_dofs(); }
        SizeType n_local_dofs() const override { return unconstrained_space()->n_local_dofs(); }

        void create_vector(Vector &v) const override { unconstrained_space()->create_vector(v); }
        void create_matrix(Matrix &m) const override { unconstrained_space()->create_matrix(m); }

        void apply_constraints(Matrix &m, const Scalar diag_value = 1.0) const override {
            if (projector_) {
                Matrix m_temp;
                projector_->restrict(m, m_temp);
                m += m_temp;

                set_zero_rows(m, constrained_indices_, diag_value);
            }

            unconstrained_space()->apply_constraints(m, diag_value);
        }

        void apply_constraints(Vector &v) const override {
            if (projector_) {
                Vector v_temp;
                projector_->restrict(v, v_temp);
                v += v_temp;

                set(v, constrained_indices_, 0.);
            }

            unconstrained_space()->apply_constraints(v);
        }

        void apply_constraints(Matrix &m, Vector &v) const override {
            if (projector_) {
                Matrix m_temp;
                projector_->restrict(m, m_temp);
                m += m_temp;

                Vector v_temp;
                projector_->restrict(v, v_temp);
                v += v_temp;

                set_zero_rows(m, constrained_indices_, 1.0);
                set(v, constrained_indices_, 0.);
            }

            unconstrained_space()->apply_constraints(m, v);
        }

        void apply_zero_constraints(Vector &vec) const override { unconstrained_space()->apply_zero_constraints(vec); }

        void add_dirichlet_boundary_condition(const std::string &name,
                                              const Scalar &value,
                                              const int component = 0) override {
            unconstrained_space()->add_dirichlet_boundary_condition(name, value, component);
        }

        bool empty() const override { return unconstrained_space()->empty(); }

        void displace(const Vector &displacement) override { unconstrained_space()->displace(displacement); }

        const std::string &name() const override { return unconstrained_space()->name(); }
        void initialize() override { unconstrained_space()->initialize(); }

        int n_var() const override { return unconstrained_space()->n_var(); }

    private:
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<Transfer> projector_;
        IndexArray constrained_indices_;
    };

}  // namespace utopia

#endif  // UTOPIA_NC_FUNCTION_SPACE_HPP

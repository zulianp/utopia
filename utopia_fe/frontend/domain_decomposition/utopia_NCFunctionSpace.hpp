#ifndef UTOPIA_NC_FUNCTION_SPACE_HPP
#define UTOPIA_NC_FUNCTION_SPACE_HPP

#include "utopia_fe_base.hpp"

#include "utopia_IPTransfer.hpp"
#include "utopia_Input.hpp"
#include "utopia_Temp.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Transfer.hpp"

#include "utopia_FunctionSpaceBase.hpp"
#include "utopia_fe_Core.hpp"

#include <memory>
#include <string>

namespace utopia {

    template <class FunctionSpace>
    class Traits<NCFunctionSpace<FunctionSpace>> : public Traits<FunctionSpace> {};

    template <class FunctionSpace>
    class NCFunctionSpace : public FunctionSpace {
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
        using Super = FunctionSpace;

        NCFunctionSpace(const Communicator &comm) : Super(comm) {}

        void read(Input &in) override {
            Super::read(in);

            std::string type;
            in.get("type", type);
            in.get("debug", debug_);

#ifdef UTOPIA_WITH_MOONOLITH
            if (type == "mortar") {
                FETransfer<FunctionSpace> transfer;
                in.get("mortar", transfer);
                if (!transfer.init_from_decomposition(make_ref(*this))) {
                    Utopia::Abort("NCFunctionSpace: failed to initialize mortar operator!");
                }

                auto temp = transfer.template build_transfer<IPTransfer>();

                Super::apply_constraints(*temp->I_ptr(), 0.0);

                Vector is_constrained = sum(temp->I(), 1);
                constrained_indices_.reserve(is_constrained.local_size());
                is_constrained.read([this](const SizeType &i, const Scalar &v) {
                    if (v > 0.99) {
                        constrained_indices_.push_back(i);
                    }
                });

                temp->I_ptr()->shift_diag(1.0);

                projector_ = temp;
            }
#endif  // UTOPIA_WITH_MOONOLITH
        }

        void restrict(const Matrix &in, Matrix &out) const {
            if (projector_) {
                Matrix temp;
                projector_->restrict(in, temp);
                // out += temp;
                out = temp;
            }
        }

        void restrict(const Vector &in, Vector &out) const {
            if (projector_) {
                projector_->restrict(in, out);
            }
        }

        void interpolate(const Vector &in, Vector &out) const {
            if (projector_) {
                if (in.is_alias(out)) {
                    Vector temp;
                    projector_->interpolate(in, temp);
                    // out += temp;
                    out = temp;
                } else {
                    projector_->interpolate(in, out);
                    // out += in;
                }
            } else {
                if (!in.is_alias(out)) {
                    out = in;
                }
            }
        }

        void project(const Vector &in, Vector &out) const {
            if (projector_) {
                projector_->interpolate(in, out);
            }
        }

        bool write(const Path &path, const Vector &x) override {
            if (projector_) {
                Vector out;
                interpolate(x, out);
                return Super::write(path, out);
            } else {
                return Super::write(path, x);
            }
        }

        void apply_constraints(Matrix &m, const Scalar diag_value = 1.0) const override {
            if (projector_) {
                Matrix m_temp;
                projector_->restrict(m, m_temp);

                if (debug_) {
                    utopia::rename("p", m_temp);
                    utopia::write("load_p.m", m_temp);
                }

                // m += m_temp;
                m = m_temp;

                utopia::set_zero_rows(m, constrained_indices_, diag_value);
            }

            Super::apply_constraints(m, diag_value);
        }

        void apply_constraints(Vector &v) const override {
            if (projector_) {
                Vector v_temp;
                projector_->restrict(v, v_temp);
                // v += v_temp;
                v = v_temp;

                utopia::set(v, constrained_indices_, 0.);
            }

            Super::apply_constraints(v);
        }

        void apply_constraints_update(Vector &v) const override { Super::apply_constraints(v); }

        void apply_constraints(Matrix &m, Vector &v) const override {
            if (projector_) {
                Matrix m_temp;
                projector_->restrict(m, m_temp);
                // m += m_temp;
                m = m_temp;

                Vector v_temp;
                projector_->restrict(v, v_temp);
                // v += v_temp;
                v = v_temp;

                set_zero_rows(m, constrained_indices_, 1.0);
                set(v, constrained_indices_, 0.);
            }

            Super::apply_constraints(m, v);
        }

        inline bool is_non_conforming() const override { return static_cast<bool>(projector_); }

        inline std::shared_ptr<Matrix> constraint_matrix() const override {
            if (projector_) {
                auto copy = std::make_shared<Matrix>(*projector_->I_ptr());
                copy->shift_diag(-1);
                return copy;
            } else {
                return nullptr;
            }
        }

    private:
        std::shared_ptr<IPTransfer> projector_;
        IndexArray constrained_indices_;
        bool debug_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_NC_FUNCTION_SPACE_HPP

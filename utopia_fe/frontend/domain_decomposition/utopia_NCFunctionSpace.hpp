#ifndef UTOPIA_NC_FUNCTION_SPACE_HPP
#define UTOPIA_NC_FUNCTION_SPACE_HPP

#include "utopia_IPTransfer.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Transfer.hpp"

#include "utopia_fe_Core.hpp"

#ifdef UTOPIA_WITH_MOONOLITH
#include "utopia_FETransferOptions.hpp"
#endif  // UTOPIA_WITH_MOONOLITH

#include <memory>

namespace utopia {

    template <class FunctionSpace>
    class NCFunctionSpace : public Configurable {
    public:
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Transfer_t = utopia::Transfer<Matrix_t, Vector_t>;
        using IPTransfer_t = utopia::IPTransfer<Matrix_t, Vector_t>;

        std::shared_ptr<FunctionSpace> space() const { return space_; }

        NCFunctionSpace(const std::shared_ptr<FunctionSpace> &space) : space_(space) {}
        NCFunctionSpace(const Communicator_t &comm) : space_(std::make_shared<FunctionSpace>(comm)) {}

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
                    Utopia::Abort("NCFunctionSpace: failed to initialize mortar");
                }

                projector_ = transfer.template build_transfer<IPTransfer_t>();
            }
#endif  // UTOPIA_WITH_MOONOLITH
        }

        void restrict(const Matrix_t &in, Matrix_t &out) const {
            if (projector_) {
                projector_->restrict(in, out);
            }
        }

        void restrict(const Vector_t &in, Vector_t &out) const {
            if (projector_) {
                projector_->restrict(in, out);
            }
        }

        // void interpolate(const Matrix_t &in, Matrix_t &out) const {
        //     out = in;
        //     if (projector_) {
        //         Vector_t temp;
        //         projector_->apply(in, temp);
        //         out += temp;
        //     }
        // }

        void interpolate(const Vector_t &in, Vector_t &out) const {
            out = in;
            if (projector_) {
                Vector_t temp;
                projector_->interpolate(in, temp);
                out += temp;
            }
        }

        void project(const Vector_t &in, Vector_t &out) const {
            if (projector_) {
                projector_->interpolate(in, out);
            }
        }

    private:
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<Transfer_t> projector_;
    };

    template <class FunctionSpace>
    class Traits<NCFunctionSpace<FunctionSpace>> : public Traits<FunctionSpace> {};

}  // namespace utopia

#endif  // UTOPIA_NC_FUNCTION_SPACE_HPP

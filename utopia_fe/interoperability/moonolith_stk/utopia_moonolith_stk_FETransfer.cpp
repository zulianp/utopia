#include "utopia_moonolith_stk_FETransfer.hpp"

#include "utopia_moonolith_FETransfer.hpp"

#include "utopia_moonolith_stk.hpp"

namespace utopia {

    namespace stk {

        class FETransfer::Impl {
        public:
            std::shared_ptr<moonolith::FunctionSpace> from;
            std::shared_ptr<moonolith::FunctionSpace> to;
            moonolith::FETransfer transfer;
            bool export_converted_mesh{false};
            bool remove_constrained_dofs{false};
        };

        void FETransfer::read(Input &in) {
            impl_->transfer.read(in);
            in.get("export_converted_mesh", impl_->export_converted_mesh);
            in.get("remove_constrained_dofs", impl_->remove_constrained_dofs);
        }

        void FETransfer::describe(std::ostream &os) const { impl_->transfer.describe(os); }

        bool FETransfer::init_from_decomposition(const std::shared_ptr<FunctionSpace> &from_and_to) {
            impl_->from = std::make_shared<moonolith::FunctionSpace>(from_and_to->comm());
            impl_->to = impl_->from;

            extract_trace_space(*from_and_to, *impl_->from);

            if (impl_->export_converted_mesh) {
                impl_->from->mesh().write("from_and_to.vtu");
            }

            return impl_->transfer.init(impl_->from);
        }

        bool FETransfer::init(const std::shared_ptr<FunctionSpace> &from_and_to) {
            impl_->from = std::make_shared<moonolith::FunctionSpace>(from_and_to->comm());
            impl_->to = impl_->from;

            convert_function_space(*from_and_to, *impl_->from);

            assert(from_and_to->mesh().n_local_elements() == impl_->from->mesh().n_local_elements());

            if (impl_->export_converted_mesh) {
                impl_->from->mesh().write("from_and_to.vtu");
            }

            if(!impl_->transfer.init(impl_->from)) {
                return false;
            }

            if(impl_->remove_constrained_dofs) {
                from_and_to->apply_constraints(*transfer_matrix(), 0);
            }

            return true;
        }

        bool FETransfer::init(const std::shared_ptr<FunctionSpace> &from, const std::shared_ptr<FunctionSpace> &to) {
            impl_->from = std::make_shared<moonolith::FunctionSpace>(from->comm());
            impl_->to = std::make_shared<moonolith::FunctionSpace>(to->comm());

            convert_function_space(*from, *impl_->from);
            convert_function_space(*to, *impl_->to);

            assert(from->mesh().n_local_elements() == impl_->from->mesh().n_local_elements());
            assert(to->mesh().n_local_elements() == impl_->to->mesh().n_local_elements());

            if (impl_->export_converted_mesh) {
                impl_->from->mesh().write("from.vtu");
                impl_->to->mesh().write("to.vtu");
            }

            if(!impl_->transfer.init(impl_->from, impl_->to)) {
                return false;
            }

            if(impl_->remove_constrained_dofs) {
                to->apply_constraints(*transfer_matrix(), 0);
            }

            return true;
        }

        void FETransfer::verbose(const bool val) { impl_->transfer.verbose(val); }

        void FETransfer::clear() {
            impl_->transfer.clear();
            impl_->from = nullptr;
            impl_->to = nullptr;
        }

        void FETransfer::set_options(const FETransferOptions &options) { impl_->transfer.set_options(options); }

        bool FETransfer::empty() const { return impl_->transfer.empty(); }

        bool FETransfer::apply(const Vector &from, Vector &to) const { return impl_->transfer.apply(from, to); }

        bool FETransfer::apply(const Matrix &to_matrix, Matrix &matrix_in_from_space) const {
            return impl_->transfer.apply(to_matrix, matrix_in_from_space);
        }

        Size FETransfer::size() const { return impl_->transfer.size(); }

        Size FETransfer::local_size() const { return impl_->transfer.local_size(); }

        FETransfer::Communicator &FETransfer::comm() { return impl_->transfer.comm(); }

        const FETransfer::Communicator &FETransfer::comm() const { return impl_->transfer.comm(); }

        bool FETransfer::apply_transpose(const Vector &from, Vector &to) const {
            return impl_->transfer.apply_transpose(from, to);
        }

        bool FETransfer::write(const Path &path) const { return impl_->transfer.write(path); }

        FETransfer::FETransfer() : impl_(utopia::make_unique<Impl>()) {}

        FETransfer::~FETransfer() = default;

        std::shared_ptr<FETransfer::Matrix> FETransfer::transfer_matrix() const {
            return impl_->transfer.transfer_matrix();
        }

    }  // namespace stk
}  // namespace utopia

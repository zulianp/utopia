#include "utopia_moonolith_libmesh_FETransfer.hpp"

#include "utopia_moonolith_FETransfer.hpp"

#include "utopia_moonolith_libmesh.hpp"

namespace utopia {

    namespace libmesh {

        class FETransfer::Impl {
        public:
            std::shared_ptr<moonolith::FunctionSpace> from;
            std::shared_ptr<moonolith::FunctionSpace> to;
            moonolith::FETransfer transfer;
        };

        void FETransfer::read(Input &in) { impl_->transfer.read(in); }

        void FETransfer::describe(std::ostream &os) const { impl_->transfer.describe(os); }

        bool FETransfer::init(const std::shared_ptr<FunctionSpace> &from, const std::shared_ptr<FunctionSpace> &to) {
            impl_->from = std::make_shared<moonolith::FunctionSpace>(from->comm());
            impl_->to = std::make_shared<moonolith::FunctionSpace>(to->comm());

            convert_function_space(*from, *impl_->from);
            convert_function_space(*to, *impl_->to);

            return impl_->transfer.init(impl_->from, impl_->to);
        }

        void FETransfer::clear() {
            impl_->transfer.clear();
            impl_->from = nullptr;
            impl_->to = nullptr;
        }

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

    }  // namespace libmesh
}  // namespace utopia

#ifndef UTOPIA_MOONOLITH_FE_TRANSFER_HPP
#define UTOPIA_MOONOLITH_FE_TRANSFER_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Operator.hpp"

#include "utopia_moonolith_ForwardDeclarations.hpp"
#include "utopia_moonolith_FunctionSpace.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_FETransferOptions.hpp"
#include "utopia_Transfer.hpp"

namespace utopia {

    namespace moonolith {

        class FETransfer : public Operator<Traits<FunctionSpace>::Vector>, public Configurable, public Describable {
        public:
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Vector = Traits<FunctionSpace>::Vector;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;
            using Communicator = Traits<FunctionSpace>::Communicator;

            void set_options(const FETransferOptions &options);

            void read(Input &in) override;
            void describe(std::ostream &) const override;

            bool init(const std::shared_ptr<FunctionSpace> &from, const std::shared_ptr<FunctionSpace> &to);
            bool init(const std::shared_ptr<FunctionSpace> &from_and_to);

            void clear();
            bool empty() const;
            bool apply(const Vector &from, Vector &to) const override;
            bool apply(const Matrix &to_matrix, Matrix &from_matrix, const bool reuse_matrix = false) const;
            Size size() const override;
            Size local_size() const override;

            void verbose(const bool val);

            Communicator &comm() override;
            const Communicator &comm() const override;

            bool apply_transpose(const Vector &from, Vector &to) const;
            bool write(const Path &) const;

            FETransfer();
            ~FETransfer();

            template <class TransferType>
            inline std::shared_ptr<TransferType> build_transfer() const {
                return std::make_shared<TransferType>(transfer_matrix());
            }

            std::shared_ptr<Matrix> transfer_matrix() const;

            void set_constraint_matrix_from(const std::shared_ptr<Matrix> &c);
            void set_constraint_matrix_to(const std::shared_ptr<Matrix> &c);

            // void rebalance_the_from_system();

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace moonolith

    template <>
    class FETransfer<utopia::moonolith::FunctionSpace> : public utopia::moonolith::FETransfer {};
}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_FE_TRANSFER_HPP
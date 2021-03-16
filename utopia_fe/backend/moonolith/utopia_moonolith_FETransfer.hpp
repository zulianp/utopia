#ifndef UTOPIA_MOONOLITH_FE_TRANSFER_HPP
#define UTOPIA_MOONOLITH_FE_TRANSFER_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Operator.hpp"

#include "utopia_moonolith_ForwardDeclarations.hpp"
#include "utopia_moonolith_FunctionSpace.hpp"

namespace utopia {

    namespace moonolith {

        class FETransfer : public Operator<Traits<FunctionSpace>::Vector>, public Configurable, public Describable {
        public:
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Vector = Traits<FunctionSpace>::Vector;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;
            using Communicator = Traits<FunctionSpace>::Communicator;

            void read(Input &in) override;
            void describe(std::ostream &) const override;

            bool init(const std::shared_ptr<FunctionSpace> &from, const std::shared_ptr<FunctionSpace> &to);
            void clear();
            bool empty() const;
            bool apply(const Vector &from, Vector &to) const override;
            Size size() const override;
            Size local_size() const override;

            Communicator &comm() override;
            const Communicator &comm() const override;

            bool apply_transpose(const Vector &from, Vector &to) const;
            bool write(const Path &) const;

            FETransfer();
            ~FETransfer();

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace moonolith
}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_FE_TRANSFER_HPP
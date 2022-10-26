#ifndef UTOPIA_MOONOLITH_STK_CONTACT_HPP
#define UTOPIA_MOONOLITH_STK_CONTACT_HPP

#include "utopia_fe_Core.hpp"
#include "utopia_ui.hpp"

#include "utopia_moonolith_Contact.hpp"
#include "utopia_stk_FunctionSpace.hpp"

namespace utopia {

    namespace stk {

        class Contact : public ContactInterface<stk::FunctionSpace> {
        public:
            using Params = utopia::moonolith::Contact::Params;
            using Mesh = utopia::stk::Mesh;
            using FunctionSpace = utopia::stk::FunctionSpace;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Vector = Traits<FunctionSpace>::Vector;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;

            void set_params(const Params &params);

            void read(Input &in) override;
            void describe(std::ostream &os) const override;
            bool assemble(FunctionSpace &space) override;
            void transform(const Matrix &in, Matrix &out) override;
            void transform(const Vector &in, Vector &out) override;
            void inverse_transform(const Vector &in, Vector &out) override;

            std::shared_ptr<Matrix> orthogonal_transformation() override;
            std::shared_ptr<Matrix> mass_matrix() override;

            Contact();
            virtual ~Contact();

            const Vector &gap() const override;
            const Vector &is_contact() const override;
            const Vector &normals() const override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace stk

    template <>
    class Contact<utopia::stk::FunctionSpace> final : public utopia::stk::Contact {};
}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_STK_CONTACT_HPP

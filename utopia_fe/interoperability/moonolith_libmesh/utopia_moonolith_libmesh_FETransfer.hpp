#ifndef UTOPIA_MOONOLITH_LIBMESH_FETRANSFER_HPP
#define UTOPIA_MOONOLITH_LIBMESH_FETRANSFER_HPP

#ifndef UTOPIA_MOONOLITH_STK_FE_TRANSFER_HPP
#define UTOPIA_MOONOLITH_STK_FE_TRANSFER_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Operator.hpp"

#include "utopia_libmesh_ForwardDeclarations.hpp"
#include "utopia_libmesh_FunctionSpace_new.hpp"

#include "utopia_fe_Core.hpp"

namespace utopia {

    namespace libmesh {

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
            bool apply(const Matrix &to_matrix, Matrix &matrix_in_from_space) const;
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

    }  // namespace libmesh

    template <>
    class FETransfer<utopia::libmesh::FunctionSpace> : public utopia::libmesh::FETransfer {};

}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_STK_FE_TRANSFER_HPP

#endif  // UTOPIA_MOONOLITH_LIBMESH_FETRANSFER_HPP
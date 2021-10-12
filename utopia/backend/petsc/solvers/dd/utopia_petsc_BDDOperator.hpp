#ifndef UTOPIA_PETSC_BDDOPERATOR_HPP
#define UTOPIA_PETSC_BDDOPERATOR_HPP

#include "utopia_Operator.hpp"

#include <memory>

namespace utopia {

    template <class Matrix, class Vector>
    class BDDOperator : public Operator<Vector> {
    public:
        using Communicator = typename Traits<Vector>::Communicator;

        BDDOperator();
        ~BDDOperator();

        bool update(const std::shared_ptr<Matrix> &matrix);

        bool apply(const Vector &x_G, Vector &rhs_G) const override;
        Size size() const override;

        Size local_size() const override;

        Communicator &comm() override;
        const Communicator &comm() const override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_BDDOPERATOR_HPP

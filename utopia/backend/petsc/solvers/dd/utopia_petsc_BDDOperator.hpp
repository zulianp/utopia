#ifndef UTOPIA_PETSC_BDDOPERATOR_HPP
#define UTOPIA_PETSC_BDDOPERATOR_HPP

#include "utopia_Input.hpp"
#include "utopia_Operator.hpp"
#include "utopia_Preconditioner.hpp"

#include <memory>
#include <vector>

namespace utopia {

    template <class Matrix, class Vector>
    class BDDOperator : public Operator<Vector>, public Configurable {
    public:
        using Communicator = typename Traits<Vector>::Communicator;
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using IndexArray = typename Traits::IndexArray;
        using IndexSet = typename Traits::IndexSet;
        using Layout = typename Traits::Layout;
        using Selector = std::vector<bool>;

        BDDOperator();
        ~BDDOperator();

        void read(Input &in) override;

        bool initialize(const std::shared_ptr<const Vector> &rhs);
        bool initialize(const std::shared_ptr<const Matrix> &matrix);
        bool finalize(const Vector &x_G, Vector &x);

        void create_vector(Vector &x_G);

        const Vector &righthand_side() const;
        void diag(Vector &d) const;

        std::shared_ptr<Preconditioner<Vector>> create_preconditioner() const;
        std::shared_ptr<Matrix> reduced_matrix() const;
        const IndexSet &skeleton_dofs() const;

        void select(const Vector &x, Vector &x_G) const;
        bool apply(const Vector &x_G, Vector &rhs_G) const override;
        Size size() const override;

        Layout vector_layout() const;

        Size local_size() const override;

        Communicator &comm() override;
        const Communicator &comm() const override;

        Selector &selector();

        UTOPIA_NVCC_PRIVATE
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_BDDOPERATOR_HPP

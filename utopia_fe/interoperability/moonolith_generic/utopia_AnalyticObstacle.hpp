#ifndef UTOPIA_ANALYTIC_OBSTACLE_HPP
#define UTOPIA_ANALYTIC_OBSTACLE_HPP

#include "utopia_ContactInterface.hpp"
#include "utopia_Field.hpp"
#include "utopia_GradientField.hpp"
#include "utopia_Transfer.hpp"

#include <memory>

namespace utopia {

    template <class FunctionSpace>
    class AnalyticObstacle : public ContactInterface<FunctionSpace> {
    public:
        using Super = utopia::ContactInterface<FunctionSpace>;
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using Comm = typename Traits<FunctionSpace>::Communicator;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using Field = utopia::Field<FunctionSpace>;
        using GradientField = utopia::GradientField<FunctionSpace>;

        void read(Input &in) override;
        void describe(std::ostream &os) const override;

        bool init(const std::shared_ptr<FunctionSpace> &domain);
        void set_transfer(const std::shared_ptr<Transfer> &transfer);

        bool assemble(FunctionSpace &space) override;
        void transform(const Matrix &in, Matrix &out) override;
        void transform(const Vector &in, Vector &out) override;
        void inverse_transform(const Vector &in, Vector &out) override;

        std::shared_ptr<Matrix> orthogonal_transformation() override;

        const Vector &gap() const override;
        const Vector &is_contact() const override;
        const Vector &normals() const override;

        AnalyticObstacle();
        ~AnalyticObstacle();

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_ANALYTIC_OBSTACLE_HPP

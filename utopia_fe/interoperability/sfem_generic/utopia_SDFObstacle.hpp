#ifndef UTOPIA_SDF_OBSTACLE_HPP
#define UTOPIA_SDF_OBSTACLE_HPP

#include "utopia_ContactInterface.hpp"
#include "utopia_Field.hpp"
#include "utopia_GradientField.hpp"
#include "utopia_Transfer.hpp"

#include <memory>

namespace utopia {

    template <class FunctionSpace>
    class SDFObstacle : public ContactInterface<FunctionSpace> {
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

        bool assemble(FunctionSpace &space) override;
        void transform(const Matrix &in, Matrix &out) override;
        void transform(const Vector &in, Vector &out) override;
        void inverse_transform(const Vector &in, Vector &out) override;

        std::shared_ptr<Matrix> orthogonal_transformation() override;

        const Vector &gap() const override;
        const Vector &is_contact() const override;
        const Vector &normals() const override;

        SDFObstacle();
        ~SDFObstacle();

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_SDF_OBSTACLE_HPP

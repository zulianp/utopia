#ifndef UTOPIA_I_OBSTACLE_HPP
#define UTOPIA_I_OBSTACLE_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class FunctionSpace>
    class IObstacle : public Configurable, public Describable {
    public:
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using Comm = typename Traits<FunctionSpace>::Communicator;

        virtual bool assemble(FunctionSpace &space) = 0;
        virtual void transform(const Matrix &in, Matrix &out) = 0;
        virtual void transform(const Vector &in, Vector &out) = 0;
        virtual void inverse_transform(const Vector &in, Vector &out) = 0;

        virtual std::shared_ptr<Matrix> orthogonal_transformation() const { return nullptr; }

        virtual const Vector &gap() const = 0;
        virtual const Vector &is_contact() const = 0;
        virtual const Vector &normals() const = 0;

        virtual ~IObstacle() = default;
    };

}  // namespace utopia

#endif  // UTOPIA_I_OBSTACLE_HPP

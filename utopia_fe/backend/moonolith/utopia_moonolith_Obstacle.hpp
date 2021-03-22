#ifndef UTOPIA_MOONOLITH_OBSTACLE_HPP
#define UTOPIA_MOONOLITH_OBSTACLE_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_moonolith_ForwardDeclarations.hpp"
#include "utopia_moonolith_FunctionSpace.hpp"

#include <memory>

namespace utopia {
    namespace moonolith {

        class Obstacle final : public Configurable, public Describable {
        public:
            using Vector = Traits<FunctionSpace>::Vector;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using Comm = Traits<FunctionSpace>::Communicator;

            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            bool init_obstacle(const Mesh &obstacle_mesh);
            bool assemble(const FunctionSpace &space);
            void transform(const Matrix &in, Matrix &out);
            void transform(const Vector &in, Vector &out);
            void inverse_transform(const Vector &in, Vector &out);

            const Vector &gap() const;
            const Vector &is_contact() const;

            Obstacle();
            ~Obstacle();

        private:
            class Impl;
            class Output;
            class Params;

            template <int Dim>
            class ImplD;

            Output &output();
            const Output &output() const;

            std::unique_ptr<Params> params_;
            std::unique_ptr<Output> output_;

            std::unique_ptr<Impl> impl_;
        };
    }  // namespace moonolith
}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_OBSTACLE_HPP
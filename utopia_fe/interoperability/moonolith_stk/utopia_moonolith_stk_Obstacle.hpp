#ifndef UTOPIA_MOONOLITH_STK_OBSTACLE_HPP
#define UTOPIA_MOONOLITH_STK_OBSTACLE_HPP

#include "utopia_fe_Core.hpp"
#include "utopia_ui.hpp"

#include "utopia_moonolith_Obstacle.hpp"
#include "utopia_stk_FunctionSpace.hpp"

namespace utopia {

    namespace stk {

        class Obstacle {
        public:
            using Params = utopia::moonolith::Obstacle::Params;
            using Mesh = utopia::stk::Mesh;
            using FunctionSpace = utopia::stk::FunctionSpace;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Vector = Traits<FunctionSpace>::Vector;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;

            void set_params(const Params &params);
            bool assemble(const FunctionSpace &space);
            bool init_obstacle(const Mesh &mesh);

            void transform(const Matrix &in, Matrix &out);
            void transform(const Vector &in, Vector &out);
            void inverse_transform(const Vector &in, Vector &out);

            Obstacle();
            virtual ~Obstacle();

            const Vector &gap() const;
            const Vector &is_contact() const;
            const Vector &normals() const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace stk

    template <>
    class Obstacle<utopia::stk::FunctionSpace> final : public utopia::stk::Obstacle {};
}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_STK_OBSTACLE_HPP

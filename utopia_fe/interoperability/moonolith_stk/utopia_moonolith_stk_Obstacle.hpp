#ifndef UTOPIA_MOONOLITH_STK_OBSTACLE_HPP
#define UTOPIA_MOONOLITH_STK_OBSTACLE_HPP

#include "utopia_ContactInterface.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_ui.hpp"

#include "utopia_moonolith_Obstacle.hpp"
#include "utopia_stk_FunctionSpace.hpp"

namespace utopia {

    namespace stk {

        class Obstacle : public ContactInterface<stk::FunctionSpace> {
        public:
            using Params = utopia::moonolith::Obstacle::Params;
            using Mesh = utopia::stk::Mesh;
            using FunctionSpace = utopia::stk::FunctionSpace;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Vector = Traits<FunctionSpace>::Vector;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;

            void set_params(const Params &params);
            bool init_obstacle(const Mesh &mesh);

            void read(Input &in) override;
            void describe(std::ostream &os) const override;
            bool assemble(FunctionSpace &space) override;
            void transform(const Matrix &in, Matrix &out) override;
            void transform(const Vector &in, Vector &out) override;
            void inverse_transform(const Vector &in, Vector &out) override;

            std::shared_ptr<Matrix> orthogonal_transformation() override;
            std::shared_ptr<Matrix> mass_matrix() override;

            Obstacle();
            virtual ~Obstacle();

            const Vector &gap() const override;
            const Vector &is_contact() const override;
            const Vector &normals() const override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace stk

    template <>
    class Obstacle<utopia::stk::FunctionSpace> final : public utopia::stk::Obstacle {};
}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_STK_OBSTACLE_HPP

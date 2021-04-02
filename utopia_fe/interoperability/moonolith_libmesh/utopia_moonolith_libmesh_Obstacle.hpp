#ifndef UTOPIA_OBSTACLE_HPP
#define UTOPIA_OBSTACLE_HPP

#include "utopia_fe_Core.hpp"
#include "utopia_ui.hpp"

#include "utopia_libmesh_FunctionSpace_new.hpp"
#include "utopia_moonolith_Obstacle.hpp"

namespace utopia {

    namespace libmesh {

        class Obstacle {
        public:
            using Params = utopia::moonolith::Obstacle::Params;
            using Mesh = utopia::libmesh::Mesh;
            using FunctionSpace = utopia::libmesh::FunctionSpace;
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

    }  // namespace libmesh

    template <>
    class Obstacle<utopia::libmesh::FunctionSpace> final : public utopia::libmesh::Obstacle {};
}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_HPP

#ifndef UTOPIA_MOONOLITH_OBSTACLE_HPP
#define UTOPIA_MOONOLITH_OBSTACLE_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_moonolith_ForwardDeclarations.hpp"
#include "utopia_moonolith_FunctionSpace.hpp"

#include <memory>
#include <unordered_set>

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
            const Vector &normals() const;

            Obstacle();
            ~Obstacle();

            // class SurfaceDescriptor {
            // public:
            //     std::string name;
            //     int id{-1};
            // };

            class Params : public Configurable {
            public:
                int variable_number{0};
                Scalar gap_negative_bound{-0.0001};
                Scalar gap_positive_bound{0.1};
                std::unordered_set<int> tags;
                bool invert_face_orientation{false};

                void read(Input &in) override;
            };

            void set_params(const Params &params);
            inline const Params &params() const { return *params_; }

        private:
            class Impl;
            class Output;

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

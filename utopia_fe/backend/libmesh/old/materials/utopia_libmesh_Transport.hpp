#ifndef UTOPIA_LIBMESH_TRANSPORT_HPP
#define UTOPIA_LIBMESH_TRANSPORT_HPP

#include "utopia_Field.hpp"

#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"
#include "utopia_libmesh_ForwardDeclarations.hpp"
#include "utopia_libmesh_FunctionSpace_new.hpp"

// Legacy
#include "utopia_Model.hpp"

namespace utopia {

    template <>
    class FEAssembler<libmesh::FunctionSpace>
        : public Model<Traits<libmesh::FunctionSpace>::Matrix, Traits<libmesh::FunctionSpace>::Vector> {
    public:
        using Matrix = Traits<libmesh::FunctionSpace>::Matrix;
        using Vector = Traits<libmesh::FunctionSpace>::Vector;

        virtual ~FEAssembler() = default;
        virtual bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) = 0;
        virtual bool assemble(const Vector &x, Vector &gradient) = 0;
        virtual bool assemble(const Vector &x, Matrix &hessian) = 0;

        virtual bool assemble(Matrix &) { return false; }

        inline bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override {
            return this->assemble(x, hessian, gradient);
        }

        inline bool assemble_hessian(const Vector &x, Matrix &hessian) override { return this->assemble(x, hessian); }
        inline bool assemble_hessian(Matrix &hessian) override { return this->assemble(hessian); }

        inline bool assemble_gradient(const Vector &x, Vector &gradient) override {
            return this->assemble(x, gradient);
        }

        inline void set_environment(const std::shared_ptr<Environment<libmesh::FunctionSpace>> &env) { env_ = env; }
        const std::shared_ptr<Environment<libmesh::FunctionSpace>> &environment() const { return env_; }

        inline void set_space(const std::shared_ptr<libmesh::FunctionSpace> &space) { space_ = space; }

        inline std::shared_ptr<libmesh::FunctionSpace> space() { return space_; }

        inline void read(Input &in) override {
            if (env_) {
                std::string space_name;
                in.get("space", space_name);

                if (!space_name.empty()) {
                    assert(!space_);

                    space_ = env_->find_space(space_name);
                }
            }
        }

        inline void clear() override {}
        virtual std::string name() const = 0;

    private:
        std::shared_ptr<Environment<libmesh::FunctionSpace>> env_;
        std::shared_ptr<libmesh::FunctionSpace> space_;
    };

    namespace libmesh {

        class Transport final : public FEAssembler<libmesh::FunctionSpace> {
        public:
            using Matrix = Traits<libmesh::FunctionSpace>::Matrix;
            using Vector = Traits<libmesh::FunctionSpace>::Vector;
            using Scalar = typename Traits<Vector>::Scalar;
            using Field = utopia::Field<libmesh::FunctionSpace>;
            using Super = utopia::FEAssembler<libmesh::FunctionSpace>;

            bool is_linear() const override { return true; }

            bool assemble(const Vector &x, Matrix &jacobian, Vector &fun) override;
            bool assemble(const Vector &x, Vector &gradient) override {
                assert(false);
                return false;
            }
            bool assemble(const Vector &x, Matrix &hessian) override {
                assert(false);
                return false;
            }
            bool assemble(Matrix &hessian) override {
                assert(false);
                return false;
            }

            void clear() override;
            void read(Input &in) override;

            inline std::string name() const override { return "Transport"; }

            Transport();
            ~Transport();

            void set_pressure_field(const std::shared_ptr<Field> &field);

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;

            void init();
            bool valid() const;
        };

        class Mass final : public FEAssembler<libmesh::FunctionSpace> {
        public:
            using Matrix = Traits<libmesh::FunctionSpace>::Matrix;
            using Vector = Traits<libmesh::FunctionSpace>::Vector;
            using Scalar = typename Traits<Vector>::Scalar;
            using Field = utopia::Field<libmesh::FunctionSpace>;
            using Super = utopia::FEAssembler<libmesh::FunctionSpace>;

            inline bool is_linear() const override { return true; }

            bool assemble(const Vector &x, Matrix &jacobian, Vector &fun) override;

            bool assemble(const Vector &x, Vector &gradient) override;
            bool assemble(const Vector &x, Matrix &hessian) override;
            bool assemble(Matrix &hessian) override;
            void clear() override;
            void read(Input &in) override;

            inline std::string name() const override { return "Mass"; }

            Mass();
            ~Mass();

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;

            void init();
            bool valid() const;
        };
    }  // namespace libmesh
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_TRANSPORT_HPP

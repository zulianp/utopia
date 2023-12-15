#ifndef UTOPIA_SOLVER_FUNCTION_HPP
#define UTOPIA_SOLVER_FUNCTION_HPP

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"

namespace utopia {
    /**
     * @brief      Base class for Nonlinear Function. All application context needed
     * by solver is usually provided inside of this functions. In optimization
     * settings, user needs to supply value(energy), gradient/
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Vector>
    class FunctionBase : public Configurable {
    public:
        DEF_UTOPIA_SCALAR(Vector);

        ~FunctionBase() override = default;

        void read(Input & /*in*/) override {}

        virtual bool value(const Vector & /*x*/, Scalar & /*value*/) const = 0;
        virtual bool gradient(const Vector & /*x*/, Vector & /*g*/) const = 0;
        virtual bool update(const Vector & /*x*/) { return true; }
        virtual bool project_onto_feasibile_region(Vector & /*x*/) const { return true; }

        virtual bool value_and_gradient(const Vector &x, Scalar &val, Vector &grad) const {
            return this->value(x, val) && this->gradient(x, grad);
        }

        virtual void create_vector(Vector &) const {}
    };

    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class Function : public FunctionBase<Vector> {
    public:
        DEF_UTOPIA_SCALAR(Vector);

        ~Function() override = default;

        virtual bool is_hessian_constant() const { return false; }

        virtual bool hessian(const Vector &x, Matrix &H) const = 0;
        virtual bool hessian(const Vector & /*point*/, Matrix & /*result*/, Matrix & /*preconditioner*/) const {
            return false;
        }

        /// If both Hessian and gradient are required at the same time use this method (it can be used to be more
        /// efficient)
        virtual bool hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const {
            return this->gradient(x, g) && this->hessian(x, H);
        }

        virtual bool hessian_and_gradient(const Vector &x, Matrix &H, Matrix &preconditioner, Vector &g) const {
            return this->gradient(x, g) && this->hessian(x, H, preconditioner);
        }

        virtual bool has_preconditioner() const { return false; }

        virtual bool initialize_hessian(Matrix & /*H*/, Matrix & /*H_pre*/) const { return false; }

        /**
         * @brief Allows the solvers to reuse allocated vectors and matrices
         */
        class Data {
        public:
            std::shared_ptr<Matrix> H;
            std::shared_ptr<Matrix> H_pre;
            std::shared_ptr<Vector> g;

            void init() {
                if (!H) {
                    H = std::make_shared<Matrix>();
                }
                if (!H_pre) {
                    H_pre = std::make_shared<Matrix>();
                }
                if (!g) {
                    g = std::make_shared<Vector>();
                }
            }
        };

        inline std::shared_ptr<Data> data() const { return data_; }

        Function() : data_(std::make_shared<Data>()) {}

    private:
        std::shared_ptr<Data> data_;
    };
}  // namespace utopia
#endif  // UTOPIA_SOLVER_FUNCTION_HPP

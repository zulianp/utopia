#ifndef UTOPIA_MARS_LAPLACE_OPERATOR_HPP
#define UTOPIA_MARS_LAPLACE_OPERATOR_HPP

#include "utopia_mars_ConcreteFEAssembler.hpp"

namespace utopia {
    namespace mars {

        template <class MarsMeshType, typename... Args>
        class LaplaceOperator : public ConcreteFEAssembler<MarsMeshType, Args...> {
        public:
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;
            using Scalar = Traits<mars::FunctionSpace>::Scalar;

            bool assemble(const Vector &x, Matrix &mat, Vector &vec) override {
                assert(false);
                return false;
            }
            bool assemble(const Vector &x, Matrix &mat) override {
                assert(false);
                return false;
            }
            bool assemble(const Vector &x, Vector &vec) override {
                assert(false);
                return false;
            }

            bool assemble(Matrix &mat) override {
                assert(false);
                return false;
            }

            void read(Input &in) override { in.get("coeff", coeff_); }

            void set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &) override {
                // assert(false);
            }

            inline bool is_linear() const override { return true; }

            Scalar coeff_{1.0};
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_LAPLACE_OPERATOR_HPP

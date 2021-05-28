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
            void read(Input &in) override { assert(false); }

            void set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &env) override {
                assert(false);
            }

            bool is_linear() const override {
                assert(false);
                return false;
            }
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_LAPLACE_OPERATOR_HPP

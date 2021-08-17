#ifndef UTOPIA_MARS_UNIFORM_FE_HPP
#define UTOPIA_MARS_UNIFORM_FE_HPP

#include "mars.hpp"

namespace utopia {
    namespace mars {

        template <class Scalar, class ExecutionSpace>
        class UniformFE {
        public:
            // template <typename Quadrature>
            // void init(Quadrature &quad) {}

        private:
            // Make texture since it is the same for every element
            // ViewMatrixType<Scalar> fun_;
            // ViewMatrixType<Scalar> grad_;
            // inv_J
            ViewMatrixType<Scalar> measure_;
        };
    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_UNIFORM_FE_HPP

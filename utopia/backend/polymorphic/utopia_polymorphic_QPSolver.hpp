#ifndef UTOPIA_POLYMORPHIC_QP_SOLVER_HPP
#define UTOPIA_POLYMORPHIC_QP_SOLVER_HPP

#include "utopia_QPSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

#include <memory>

namespace utopia {

    template<class Matrix, class Vector>
    class PolymorphicQPSolver : public QPSolver<Matrix, Vector>  {
    public:
            typedef UTOPIA_SCALAR(Vector)    			 Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector) 			 SizeType;
            typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
            typedef utopia::QPSolver<Matrix, Vector> 	 Super;
            typedef utopia::BoxConstraints<Vector>       BoxConstraints;

        public:

            PolymorphicQPSolver();
            ~PolymorphicQPSolver();
            PolymorphicQPSolver * clone() const override;
            bool apply(const Vector &rhs, Vector &sol) override;
            void read(Input &in) override;

        private:
            std::unique_ptr<QPSolver<Matrix, Vector>> impl_;

    };

}

#endif //UTOPIA_POLYMORPHIC_QP_SOLVER_HPP

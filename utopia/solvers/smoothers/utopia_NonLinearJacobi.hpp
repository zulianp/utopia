#ifndef UTOPIA_NONLINEAR_JACOBI_SMOOTHER_HPP
#define UTOPIA_NONLINEAR_JACOBI_SMOOTHER_HPP

#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"



namespace utopia
{

    /**
     * @brief      Nonlinear Jacobi smoother. Used in nonlinear MG as well as in FAS.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class NonLinearJacobi final : public NonLinearSmoother<Matrix, Vector>
    {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        typedef utopia::NonLinearSmoother<Matrix, Vector> Smoother;
        typedef utopia::Function<Matrix, Vector> Function;

    public:
        NonLinearJacobi() = default;

        bool smooth(Function & fun,  Vector &x, const Vector &rhs) override
        {
            Vector g = local_zeros(local_size(x));

            for(int i =0; i < this->sweeps(); i++)
            {
                Matrix H;
                fun.hessian(x, H);
                fun.gradient(x, g);
                g -= rhs;

                Vector d = 1./diag(H);
                H = diag(d);

                x = x - (this->relaxation_parameter() * H * g);
            }

            return true;

        }
    };

}

#endif //UTOPIA_NONLINEAR_JACOBI_SMOOTHER_HPP


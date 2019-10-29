#include "utopia_NewNeoHookean.hpp"

#include "utopia_FEEval_Local.hpp"
#include "utopia_FEEval_MultiTensor.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_AssembleLocal.hpp"

namespace utopia {

    template<class FunctionSpace, class Matrix, class Vector>
    bool NewNeoHookean<FunctionSpace, Matrix, Vector>::assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient)
    {
        using FE = utopia::FiniteElement<FunctionSpace>;

        // auto mu     = params_.var_mu();
        // auto lambda = params_.var_lambda();

        auto mu     = params_.var_mu().default_value_;
        auto lambda = params_.var_lambda().default_value_;

        //linear
        FormMatrixd gu,  alpha_gu; 
        MultiScalard dx;

        //non-linear
        MultiMatrixd guk, F, F_inv, F_inv_t, P;
        MultiScalard J, log_J, lambda_log_J, beta;
        FormMatrixd FG, FGF, stress_lin;
        FormScalard inner_Fg, PdotG;

        assemble(
            V_,
            inner(grad(trial(V_)), grad(trial(V_))) * dX,
            hessian,
            [&](FE &element, USerialMatrix &mat)
            {
                //Symbolic
                auto u = trial(element);
                
                //Symbolic to Numeric
                gu = grad(u);
                guk = grad(interpolate(x, u));
                dx = measure(element);

                //Numeric
                F = identity() + guk;
                F_inv = inv(F);
                F_inv_t = transpose(F_inv);

                J = det(F);
                log_J = logn(J);
                lambda_log_J = (rescaling_ * lambda) * log_J;

                //////////////////////////////////////////
                //Stress linearization

                FG = F_inv_t * gu;
                FGF = FG * F_inv_t;

                alpha_gu = (rescaling_ * mu) * gu;

                beta = (rescaling_ * (lambda  * log_J - mu));
                inner_Fg = inner(F_inv_t, gu);
                stress_lin = alpha_gu - beta * FGF + lambda * (inner_Fg * F_inv_t);
                
                //////////////////////////////////////////
                //Form assembly
                mat = form2(stress_lin, gu, dx);
            }
        );

        return stress(x, gradient);
    }
    
    template<class FunctionSpace, class Matrix, class Vector>
    bool NewNeoHookean<FunctionSpace, Matrix, Vector>::stress(const Vector &x, Vector &result)
    {
        using FE = utopia::FiniteElement<FunctionSpace>;

        // auto mu     = params_.var_mu();
        // auto lambda = params_.var_lambda();default_value_
        auto mu     = params_.var_mu().default_value_;
        auto lambda = params_.var_lambda().default_value_;

        //Buffers
        //linear
        FormMatrixd gu, alpha_gu; 
        MultiScalard dx;

        //non-linear
        MultiMatrixd guk, F, F_inv, F_inv_t, P;
        MultiScalard J, log_J, lambda_log_J, beta;
        FormMatrixd FG, FGF, stress_lin;
        FormScalard inner_Fg, PdotG;

        assemble(
            V_,
            grad(trial(V_)) * dX,
            result,
            [&](FE &element, USerialVector &vec)
            {
                //Symbolic
                auto u = trial(element);
                
                //Symbolic to Numeric
                gu = grad(u);
                guk = grad(interpolate(x, u));
                dx = measure(element);

                //Numeric
                F = identity() + guk;
                F_inv = inv(F);
                F_inv_t = transpose(F_inv);

                J = det(F);
                log_J = logn(J);
                lambda_log_J = (rescaling_ * lambda) * log_J;
                
                //////////////////////////////////////////
                //First Piola-Kirchoff tensor

                //axpy with scalar a
                P = F - F_inv_t;
                P *= (rescaling_ * mu);

                //axpy with multi-scalar a
                P += lambda_log_J * F_inv_t;
                PdotG = inner(P, gu);
                //////////////////////////////////////////
                //Form assembly
                vec = form1(PdotG, dx);
            }
        );

        return true;
    }

    template class NewNeoHookean<ProductFunctionSpace<LibMeshFunctionSpace>, USparseMatrix, UVector>;
}

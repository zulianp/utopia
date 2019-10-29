#include "utopia_NewNeoHookean.hpp"

#include "utopia_FEEval_Local.hpp"
#include "utopia_FEEval_MultiTensor.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_AssembleLocal.hpp"

namespace utopia {

    template<class Fun>
    inline void loop(const SizeType n, Fun f)
    {
        for(SizeType i = 0; i < n; ++i) {
            f(i);
        }
    }

    template<class FunctionSpace, class Matrix, class Vector>
    bool NewNeoHookean<FunctionSpace, Matrix, Vector>::assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient)
    {
        using FE = utopia::FiniteElement<FunctionSpace>;

        //FIXME
        // auto mu     = params_.var_mu();
        // auto lambda = params_.var_lambda();

        //x element
        MultiScalard dx;
        FormMatrixd gu;
        MultiMatrixd guk;

        //x quad-point
        USerialMatrix F, F_inv, F_inv_t, P;
        Scalar J, log_J, lambda_log_J, beta;
        USerialMatrix alpha_gu, FG, FGF, stress_lin;
        Scalar inner_Fg, PdotG;

        SizeType n_qp_max = 0;

        assemble(
            V_,
            inner(grad(trial(V_)), grad(trial(V_))) * dX,
            hessian,
            [&](FE &element, USerialMatrix &mat)
            {
                //Symbolic
                auto u = trial(element);

                auto e_id = element.ctx().block_id();
                auto mu     = params_.var_mu().get(e_id);
                auto lambda = params_.var_lambda().get(e_id);
                
                //Symbolic to Numeric
                gu = grad(u);
                guk = grad(interpolate(x, u));
                dx = measure(element);

                n_qp_max = std::max(SizeType(dx.size()), n_qp_max);

                const SizeType n_funs = gu.size();

                if(mat.rows() != n_funs || mat.cols() != n_funs) {
                    mat.resize(n_funs, n_funs);
                }
                
                mat.set(0.0);
                
                loop(dx.size(), [&](const SizeType &q) {

                    F = identity() + guk[q];
                    F_inv = inv(F);
                    F_inv_t = transpose(F_inv);
                    J = det(F);

                    log_J = std::log(J);
                    lambda_log_J = (rescaling_ * lambda) * log_J;
                    beta = (rescaling_ * (lambda  * log_J - mu));

                    for(SizeType i = 0; i < n_funs; ++i) {
                        FG = F_inv_t * gu[i][q];
                        FGF = FG * F_inv_t;
                        alpha_gu = (rescaling_ * mu) * gu[i][q];
                        inner_Fg = inner(F_inv_t, gu[i][q]);
                        stress_lin = alpha_gu - beta * FGF + lambda * (inner_Fg * F_inv_t);

                        for(SizeType j = 0; j < n_funs; ++j) {
                             const Scalar v = inner(stress_lin, gu[j][q]) * dx[q];
                             mat.add(i, j, v);
                        }
                    }
                });
            }
        );

        std::cout << "n_qp_max: " << n_qp_max << std::endl;
        return stress(x, gradient);
    }
    
    template<class FunctionSpace, class Matrix, class Vector>
    bool NewNeoHookean<FunctionSpace, Matrix, Vector>::stress(const Vector &x, Vector &result)
    {
        using FE = utopia::FiniteElement<FunctionSpace>;
        
        FormMatrixd gu, alpha_gu; 
        MultiScalard dx;
        MultiMatrixd guk;

        USerialMatrix F, F_inv, F_inv_t, P;
        Scalar J, log_J, lambda_log_J, beta;
        USerialMatrix FG, FGF, stress_lin;
        Scalar inner_Fg, PdotG;

        assemble(
            V_,
            inner(grad(trial(V_)), grad(trial(V_))) * dX,
            result,
            [&](FE &element, USerialVector &vec)
            {
                //Symbolic
                auto u = trial(element);

                auto e_id = element.ctx().block_id();
                auto mu     = params_.var_mu().get(e_id);
                auto lambda = params_.var_lambda().get(e_id);
                
                //Symbolic to Numeric
                gu = grad(u);
                guk = grad(interpolate(x, u));
                dx = measure(element);

                const SizeType n_funs = gu.size();

                if(vec.size() != n_funs) { vec.resize(n_funs); }


                vec.set(0.0);
                
                loop(dx.size(), [&](const SizeType &q) {
                    F = identity() + guk[q];
                    F_inv = inv(F);
                    F_inv_t = transpose(F_inv);
                    J = det(F);
                    log_J = std::log(J);
                    lambda_log_J = (rescaling_ * lambda) * log_J;

                    for(SizeType i = 0; i < n_funs; ++i) {
                        P = F - F_inv_t;
                        P *= (rescaling_ * mu);

                        //axpy with multi-scalar a
                        P += lambda_log_J * F_inv_t;
                        const Scalar v = inner(P, gu[i][q]) * dx[q];
                        vec.add(i, v);
                    }
                });
            }
        );

        return true;
    }

    template class NewNeoHookean<ProductFunctionSpace<LibMeshFunctionSpace>, USparseMatrix, UVector>;
}

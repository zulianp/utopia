#include "utopia_NewNeoHookean.hpp"

#include "utopia_FEEval_Local.hpp"
#include "utopia_FEEval_MultiTensor.hpp"
#include "utopia_LocalMaterial.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_AssembleLocal.hpp"

namespace utopia {

    template<class FE, class Vector>
    class LocalNeoHookean final : public LocalMaterial<FE> {
    public:
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        LocalNeoHookean(const LameeParameters &params, const Vector &x, const Scalar &rescaling)
        : params(params), x(x), rescaling(rescaling)
        {}

        void init(FE &element) override
        {
            //Symbolic
            auto u = trial(element);
            auto e_id = element.ctx().block_id();

            mu     = params.var_mu().get(e_id);
            lambda = params.var_lambda().get(e_id);

             //Symbolic to Numeric
            gu  = grad(u);
            guk = grad(interpolate(x, u));
            dx  = measure(element);
        }

        inline void assemble(FE &element, USerialMatrix &mat, USerialVector &vec) override
        {
            assemble(element, mat);
            assemble(element, vec);
        }

        void assemble(FE &element, USerialMatrix &mat) override
        {
            const SizeType n_funs = gu.size();

            loop(dx.size(), [&](const SizeType &q) {

                F = identity() + guk[q];
                F_inv = inv(F);
                F_inv_t = transpose(F_inv);
                J = det(F);

                log_J = std::log(J);
                lambda_log_J = (rescaling * lambda) * log_J;
                beta = (rescaling * (lambda  * log_J - mu));

                for(SizeType i = 0; i < n_funs; ++i) {
                    FG = F_inv_t * transpose(gu[i][q]);
                    FGF = FG * F_inv_t;
                    alpha_gu = (rescaling * mu) * gu[i][q];
                    inner_Fg = inner(F_inv_t, gu[i][q]);
                    stress_lin = alpha_gu - beta * FGF + lambda * (inner_Fg * F_inv_t);

                    //exploit symmetry
                    mat.add(i, i, inner(stress_lin, gu[i][q]) * dx[q]);

                    for(SizeType j = i+1; j < n_funs; ++j) {
                        const Scalar v = inner(stress_lin, gu[j][q]) * dx[q];
                        mat.add(i, j, v);
                        mat.add(j, i, v);
                    }
                }
            });
        }

        void assemble(FE &element, USerialVector &vec) override
        {
            const SizeType n_funs = gu.size();
            
            loop(dx.size(), [&](const SizeType &q) {
                F = identity() + guk[q];
                F_inv = inv(F);
                F_inv_t = transpose(F_inv);
                J = det(F);
                log_J = std::log(J);
                lambda_log_J = (rescaling * lambda) * log_J;

                P = F - F_inv_t;
                P *= (rescaling * mu);
                P += lambda_log_J * F_inv_t;

                for(SizeType i = 0; i < n_funs; ++i) {
                    const Scalar v = inner(P, gu[i][q]) * dx[q];
                    vec.add(i, v);
                }
            });
        }

    private:
        const LameeParameters &params;
        const Vector &x;
        const Scalar rescaling;

        //FIXME set-up from FE type
        //////////// BUFFERS ////////////////
        MultiScalard dx;
        FormMatrixd gu;
        MultiMatrixd guk;

         //x quad-point
        USerialMatrix F, F_inv, F_inv_t, P;
        Scalar J, log_J, lambda_log_J, beta, mu, lambda;
        USerialMatrix alpha_gu, FG, FGF, stress_lin;
        Scalar inner_Fg, PdotG;
        /////////////////////////////////////
    };

    template<class FunctionSpace, class Matrix, class Vector>
    bool NewNeoHookean<FunctionSpace, Matrix, Vector>::assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient)
    {
        using FE = utopia::FiniteElement<FunctionSpace>;

        LocalNeoHookean<FE, UVector> assembler(params_, x, rescaling_);

        assemble(
            V_,
            inner(grad(trial(V_)), grad(trial(V_))) * dX,
            hessian,
            [&](FE &element, USerialMatrix &mat)
            {
                assembler.init(element);
                assembler.assemble(element, mat);
            }
        );

        return stress(x, gradient);
    }
    
    template<class FunctionSpace, class Matrix, class Vector>
    bool NewNeoHookean<FunctionSpace, Matrix, Vector>::stress(const Vector &x, Vector &result)
    {
        using FE = utopia::FiniteElement<FunctionSpace>;

        LocalNeoHookean<FE, UVector> assembler(params_, x, rescaling_);

        assemble(
            V_,
            inner(grad(trial(V_)), grad(trial(V_))) * dX,
            result,
            [&](FE &element, USerialVector &vec)
            {
                assembler.init(element);
                assembler.assemble(element, vec);
            }
        );
        
        return true;
    }

    template class NewNeoHookean<ProductFunctionSpace<LibMeshFunctionSpace>, USparseMatrix, UVector>;
}

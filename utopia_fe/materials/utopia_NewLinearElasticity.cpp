#include "utopia_NewLinearElasticity.hpp"

#include "utopia_UIScalarSampler.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_FEFilter.hpp"
#include "utopia_FEEval_Filter.hpp"
#include "utopia_VonMisesStress.hpp"

#include <libmesh/tensor_value.h>
#include <libmesh/fe.h>

#include "utopia_FEEval_Local.hpp"
#include "utopia_FEEval_MultiTensor.hpp"
#include "utopia_LocalMaterial.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_AssembleLocal.hpp"

namespace utopia {

    template<class FE>
    class LocalLinearElasticity final : public LocalMaterial<FE> {
    public:
        using Scalar   = typename Traits<FE>::Scalar;
        // using SizeType = typename Traits<FE>::SizeType;

        LocalLinearElasticity(const LameeParameters &params, const Scalar &rescaling)
        : params(params), rescaling(rescaling)
        {}

        void init(FE &element) override
        {
            //Symbolic
            auto u = trial(element);
            auto e_id = element.ctx().block_id();

            mu     = params.var_mu().get(e_id);
            lambda = params.var_lambda().get(e_id);

             //Symbolic to Numeric
            grad_u = grad(u);
            div_u  = div(u);
            dx     = measure(element);

            strain.resize(grad_u.size());
        }

        void assemble(FE &element, USerialMatrix &mat) override
        {
            const SizeType n_funs = grad_u.size();

            const Scalar r_mu     = (rescaling * mu) * 0.5;
            const Scalar r_lambda = (rescaling * lambda);

            loop(dx.size(), [&](const SizeType &q) {

                const Scalar r_mu_dx     = r_mu     * dx[q];
                const Scalar r_lambda_dx = r_lambda * dx[q];

                for(SizeType i = 0; i < n_funs; ++i) {
                    strain[i] = transpose(grad_u[i][q]) + grad_u[i][q];
                }

                for(SizeType i = 0; i < n_funs; ++i) {
                    const auto &div_i = div_u[i][q];

                    //exploit symmetry
                    const Scalar v_ii = 
                            r_mu_dx     * inner(strain[i] , strain[i]) +
                            r_lambda_dx * inner(div_i, div_i);

                    mat.add(i, i, v_ii);

                    for(SizeType j = i+1; j < n_funs; ++j) {
                        
                        const Scalar v_ij = 
                            r_mu_dx     * inner(strain[i] , strain[j]) +
                            r_lambda_dx * inner(div_i, div_u[j][q]);

                        mat.add(i, j, v_ij);
                        mat.add(j, i, v_ij);
                    }
                }
            });
        }

        void assemble(FE &, USerialVector &) override
        {

        }

    private:
        const LameeParameters &params;
        const Scalar rescaling;

        //FIXME set-up from FE type
        //////////// BUFFERS ////////////////
        MultiScalard dx;
        FormMatrixd grad_u;
        MultiMatrixd strain;
        FormScalard div_u;
        Scalar mu, lambda;

    };

    template<class FunctionSpaceT, class Matrix, class Vector>
    NewLinearElasticity<FunctionSpaceT, Matrix, Vector>::NewLinearElasticity(FunctionSpaceT &V, const LameeParameters &params)
    : V_(V), params_(params), initialized_(false), rescaling_(1.0)
    {}

    template<class FunctionSpaceT, class Matrix, class Vector>
    NewLinearElasticity<FunctionSpaceT, Matrix, Vector>::~NewLinearElasticity() {}

    template<class FunctionSpaceT, class Matrix, class Vector>
    bool NewLinearElasticity<FunctionSpaceT, Matrix, Vector>::normal_stress(const UVector &x, UVector &out, const int subspace)
    {
        auto u  = trial(V_);
        auto vx = test(V_[subspace]);
        
        auto mu     = params_.var_mu();
        auto lambda = params_.var_lambda();

        auto uk = interpolate(x, u);

        auto strain = transpose(grad(uk)) + grad(uk);
        auto stress = mu * strain + lambda * trace(strain) * identity();
        auto normal_stress = dot(normal(), stress * normal());
        
        UVector mass_vector;
        bool ok = utopia::assemble(surface_integral(dot(normal_stress, vx)), out); assert(ok);
        if(!ok) return false;

        utopia::assemble(surface_integral(dot(coeff(1.0), vx)), mass_vector);

        e_pseudo_inv(mass_vector, mass_vector, 1e-14);
        out = e_mul(mass_vector, out);
        return true;
    }  

    template<class FunctionSpaceT, class Matrix, class Vector>
    bool NewLinearElasticity<FunctionSpaceT, Matrix, Vector>::von_mises_stress(const UVector &x, UVector &out, const int subspace)
    {
        auto u = trial(V_);
        auto vx = test(V_[subspace]);
        
        auto mu     = params_.var_mu();
        auto lambda = params_.var_lambda();

        auto uk = interpolate(x, u);

        auto strain = transpose(grad(uk)) + grad(uk);
        auto stress = mu * strain + lambda * trace(strain) * identity();

        auto vm = filter(stress, VonMisesStress::apply);

        UVector mass_vector;
        bool ok = utopia::assemble(integral(dot(vm, vx)), out); assert(ok);
        if(!ok) return false;

        utopia::assemble(integral(dot(coeff(1.0), vx)), mass_vector);

        e_pseudo_inv(mass_vector, mass_vector, 1e-14);
        out = e_mul(mass_vector, out);
        return true;
    }

    template<class FunctionSpaceT, class Matrix, class Vector>
    bool NewLinearElasticity<FunctionSpaceT, Matrix, Vector>::assemble_hessian(Matrix &hessian)
    {
        using FE = utopia::FiniteElement<FunctionSpaceT>;
        LocalLinearElasticity<FE> assembler(params_, rescaling_);

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
        
        return true;
    }

    template class NewLinearElasticity<ProductFunctionSpace<LibMeshFunctionSpace>, USparseMatrix, UVector>;

}

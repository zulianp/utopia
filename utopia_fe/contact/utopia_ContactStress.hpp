#ifndef UTOPIA_CONTACT_STRESS_HPP
#define UTOPIA_CONTACT_STRESS_HPP

#include "utopia_LameeParameters.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_QuadratureUtils.hpp"
#include "utopia_LinearElasticity.hpp"

#include "utopia_TransferUtils.hpp"

namespace utopia {

    template<class FunctionSpaceT, class Matrix, class Vector>
    class ContactStress final {
    public:
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;
        typedef typename TraitsT::Vector ElementVector;

        ContactStress(FunctionSpaceT &V, const LameeParameters &params)
        : V_(V), params_(params), elast_(P1_, params_)
        {
            init();
        }

        void init();

        bool assemble(const UVector &x, UVector &result);

        // template<class Expr>
        // void init_context_on(const Expr &expr, const libMesh::Elem &e)
        // {
        //     auto quad = QuadratureUtils::nodal_quad_points(e);
        //     ctx_.set_current_element(e.id());
        //     ctx_.set_has_assembled(false);
        //     ctx_.init(expr, quad);
        // }

        // template<class Expr>
        // void reinit_context_on(const Expr &expr, const libMesh::Elem &e)
        // {
        //     auto quad = QuadratureUtils::nodal_quad_points(e);
        //     ctx_.set_current_element(e.id());
        //     ctx_.set_has_assembled(false);
        //     ctx_.reinit(expr, quad);
        // }

        // void assemble_old(const UVector &x, UVector &result)
        // {
        //     const auto &dof_map = V_.subspace(0).dof_map();
        //     auto &m = V_.subspace(0).mesh();
        //     auto u  = trial(V_);
        //     auto uk = interpolate(x, u);

            
        //     auto stress_form = params_.var_mu() * (grad(uk) + transpose(grad(uk))) + (params_.var_lambda() * div(uk)) * identity();

        //     // auto ret = FEBackend<LIBMESH_TAG>::div(uk, ctx_);

        //     Vector temp_vec = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());

        //     {
        //         Write<Vector> w_v(temp_vec, utopia::GLOBAL_ADD);
        //         ElementVector el_vec;

        //         if(elements_begin(m) != elements_end(m)) {
        //             init_context_on(stress_form, **elements_begin(m));

        //             for(auto it = elements_begin(m); it != elements_end(m); ++it) {
        //                 const auto &e = **it;

        //                 if(it != elements_begin(m)) {
        //                     reinit_context_on(stress_form, e);
        //                 }

        //                 el_vec.implementation().zero();

        //                 auto cauchy_stress = quad_eval(stress_form, ctx_);

        //                 std::vector<libMesh::dof_id_type> dof_indices;
        //                 // dof_map.dof_indices(*it, dof_indices);



        //                 if(ctx_.has_assembled()) {
        //                     add_vector(el_vec.implementation(), dof_indices, temp_vec);
        //                 }
        //             }
        //         }
        //     }
        // }

    private:
        FunctionSpaceT &V_;
        FunctionSpaceT P1_;
        LameeParameters params_;
        AssemblyContext<LIBMESH_TAG> ctx_;

        USparseMatrix P1toV_, VtoP1_;

        LinearElasticity<FunctionSpaceT, Matrix, Vector> elast_;
        UVector x_p1_, stress_p1_;

        UVector inverse_mass_vector_;

        UVector p1_is_normal_component_;
        UVector p1_normals_;
        USparseMatrix p1_mat_;
};

}

#endif //UTOPIA_CONTACT_STRESS_HPP

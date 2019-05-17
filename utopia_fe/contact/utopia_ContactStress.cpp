#include "utopia_ContactStress.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_libmesh_FEBackend.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_NormalTangentialCoordinateSystem.hpp"
#include "utopia.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_fe_core.hpp"


namespace utopia {

    template<class FunctionSpaceT, class Matrix, class Vector>
    bool ContactStress<FunctionSpaceT, Matrix, Vector>::assemble(const UVector &x, UVector &result)
    {
        const auto &p1_dof_map = P1_[0].dof_map();
        x_p1_ = ghosted(p1_dof_map.n_local_dofs(), p1_dof_map.n_dofs(), p1_dof_map.get_send_list());
        x_p1_ = VtoP1_ * x;

        // bool ok = elast_.stress(x_p1_, stress_p1_); assert(ok);
        // UVector sigma = p1_mat_ * e_mul(inverse_mass_vector_, stress_p1_);

        UVector sigma;
        // bool ok = elast_.normal_stress(x_p1_, sigma);
        bool ok = elast_.von_mises_stress(x_p1_, sigma);

        result = P1toV_ * sigma;
        return ok;
    }

   template<class FunctionSpaceT, class Matrix, class Vector>
   void ContactStress<FunctionSpaceT, Matrix, Vector>::init()
    {
        auto &stress_system = V_[0].equation_systems().template add_system<libMesh::LinearImplicitSystem>("p1_gradient");

        for(auto i = 0; i < V_.n_subspaces(); ++i) {
           P1_ *= LibMeshFunctionSpace(stress_system, stress_system.add_variable("grad_" + std::to_string(i), libMesh::Order(1), libMesh::LAGRANGE) );
        }

        for(auto i = 0; i < V_.n_subspaces(); ++i) {
            P1_.subspace(i).initialize();
        }

        USparseMatrix B, D;
        // bool ok = assemble_interpolation(P1_.subspace(0), V_.subspace(0), P1toV_, D, V_.n_subspaces());  assert(ok);
        bool ok = assemble_interpolation(V_.subspace(0).mesh(), P1_.subspace(0).dof_map(), V_.subspace(0).dof_map(), P1toV_); assert(ok);
             ok = assemble_projection(V_.subspace(0), P1_.subspace(0), B, D, true, V_.n_subspaces());                         assert(ok);
             // ok = assemble_interpolation(V_.subspace(0).mesh(), V_.subspace(0).dof_map(), P1_.subspace(0).dof_map(), VtoP1_);

        UVector lumped = sum(D, 1);
        inverse_mass_vector_ = 1./lumped;

        VtoP1_ = USparseMatrix(diag(inverse_mass_vector_)) * B;
        

        USparseMatrix mass_mat;
        utopia::assemble(surface_integral(inner(trial(P1_), test(P1_))), mass_mat);

        lumped = sum(mass_mat, 1);
        e_pseudo_inv(lumped, inverse_mass_vector_, 1e-10);

        const std::vector<int> bt;
        assemble_normal_tangential_transformation(P1_[0].mesh(),
                                                  P1_[0].dof_map(),
                                                  bt,
                                                  p1_is_normal_component_,
                                                  p1_normals_,
                                                  p1_mat_);
    }

    template class ContactStress<ProductFunctionSpace<LibMeshFunctionSpace>, USparseMatrix, UVector>;
}

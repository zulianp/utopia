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
        bool ok = true;

        const auto &p1_dof_map = P1_[0].dof_map();
        x_p1_ = ghosted(p1_dof_map.n_local_dofs(), p1_dof_map.n_dofs(), p1_dof_map.get_send_list());
        x_p1_ = VtoP1_ * x;

        UVector von_mises, normal_stress, sigma;
        ok = elast_.von_mises_stress(x_p1_, von_mises, 0);  assert(ok);
        ok = elast_.normal_stress(x_p1_, normal_stress, 1); assert(ok);

        sigma = von_mises + normal_stress;

        convert(sigma, *P1_[0].equation_system().solution);
        P1_[0].equation_system().solution->close();

        ok = elast_.stress(x_p1_, sigma); assert(ok);
        result = P1toV_ * e_mul(inverse_mass_vector_, sigma);
        return ok;
    }

   template<class FunctionSpaceT, class Matrix, class Vector>
   void ContactStress<FunctionSpaceT, Matrix, Vector>::init()
    {
        auto &stress_system = V_[0].equation_systems().template add_system<libMesh::LinearImplicitSystem>("p1_gradient");

        P1_ *= LibMeshFunctionSpace(stress_system, stress_system.add_variable("von_mises", libMesh::Order(1), libMesh::LAGRANGE) );
        P1_ *= LibMeshFunctionSpace(stress_system, stress_system.add_variable("normal_stress", libMesh::Order(1), libMesh::LAGRANGE) );
        
        for(auto i = 2; i < V_.n_subspaces(); ++i) {
           P1_ *= LibMeshFunctionSpace(stress_system, stress_system.add_variable("temp_" + std::to_string(i), libMesh::Order(1), libMesh::LAGRANGE) );
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
    }

    template class ContactStress<ProductFunctionSpace<LibMeshFunctionSpace>, USparseMatrix, UVector>;
}

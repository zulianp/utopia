#include "utopia_RefactoredContactTest.hpp"
#include "utopia_ContactAssembler.hpp"

#include "utopia_ui.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIContactParams.hpp"
#include "utopia_UIMaterial.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_polymorphic_QPSolver.hpp"
#include "utopia_ContactStabilizedNewmark.hpp"
#include "utopia_InputParameters.hpp"

#include <vector>
#include <memory>

namespace utopia {

    void RefactoredContactTest::run(Input &in) {

        using ProductSpaceT    = utopia::ProductFunctionSpace<LibMeshFunctionSpace>;
        using MaterialT        = utopia::UIMaterial<ProductSpaceT, USparseMatrix, UVector>;
        using ForcingFunctionT = UIForcingFunction<ProductSpaceT, UVector>;

        in.get("contact-problem", [&](Input &in) { 
            UIMesh<libMesh::DistributedMesh> mesh(this->comm());
            UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));
            UIContactParams params;
            std::shared_ptr< ElasticMaterial<USparseMatrix, UVector> > model;

            in.get("mesh", mesh);
            in.get("space", space);
            in.get("contact", params);

            model = make_unique<MaterialT>(space.space());

            auto &V = space.space().subspace(0);

            ContactAssembler assembler;
            if(assembler.assemble(
                mesh.mesh(),
                V.dof_map(),
                params.contact_params)) {

                write("gap.e", V, assembler.gap());
                //write("is_contact.e", V, contact_data.dof_wise.is_contact);
                //write("normal.e", V, contact_data.dof_wise.normal);
            }



        });
        
    }
}


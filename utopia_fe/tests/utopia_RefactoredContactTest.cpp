#include "utopia_RefactoredContactTest.hpp"


#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_Newmark.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_ContactStabilizedNewmark.hpp"
#include "utopia_ui.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIContactParams.hpp"
#include "utopia_UIMaterial.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_InputParameters.hpp"
#include "utopia_polymorphic_QPSolver.hpp"

#include "libmesh/mesh_refinement.h"

#include "utopia_ContactAssembler.hpp"

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
            // auto forcing_function = make_unique<ForcingFunctionT>(space.space());

            // in.get("model", *model);

            auto &V = space.space().subspace(0);

            LibMeshFunctionSpaceAdapter adapter;
           
            adapter.extract_surface_init(
                make_ref(V.mesh()),
                V.dof_map(),
                params.contact_params.variable_number
            );

            adapter.print_tags();

            auto cm = std::make_shared<LibMeshCollectionManagerT>(comm()); 

            using AlogrithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<3, LibMeshFunctionSpaceAdapter>;
            using Adapter    = AlogrithmT::Adapter;

            moonolith::Communicator m_comm(comm().get());
            AlogrithmT algo(m_comm, cm);

            if(!algo.init(
                adapter,
                params.contact_params.contact_pair_tags,
                params.contact_params.search_radius
            )) {
                assert(false);
                std::cerr <<  "[Error] tree empty" << std::endl;
                return;
            }

            algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
                // auto v = isect.compute(master, slave);

                // if(v == 0.) return false;

                // vol += v;
                return true;
            });



        });
        
    }
}


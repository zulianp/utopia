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


    template<int Dim>
    bool run_contact(const ContactParams &params, LibMeshFunctionSpaceAdapter &adapter)
    {
        using AlogrithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, LibMeshFunctionSpaceAdapter>;
        using Adapter    = typename AlogrithmT::Adapter;

        auto cm = std::make_shared<LibMeshCollectionManagerT>(adapter.comm()); 

        moonolith::Communicator m_comm(cm->comm.get());

        moonolith::SearchSettings s;
        s.verbosity_level = 3;
        AlogrithmT algo(m_comm, cm, s);

        algo.init(
            adapter,
            params.contact_pair_tags,
            params.search_radius
        );

        bool ok = false;
        algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
            auto &m_space = master.collection();
            auto &m_elem  = master.elem();

            auto &s_space = slave.collection();
            auto &s_elem  = slave.elem();

            auto m_id = m_elem.id();
            auto s_id = s_elem.id();

            std::cout << m_id << "(" << master.tag() << ") -> " << s_id  << "(" << slave.tag() << ")" << std::endl;


            // auto v = isect.compute(master, slave);

            // if(v == 0.) return false;

            // vol += v;
            ok = true;
            return true;
        });

        return ok;
    }

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

            LibMeshFunctionSpaceAdapter adapter;
           
            adapter.extract_surface_init(
                make_ref(V.mesh()),
                V.dof_map(),
                params.contact_params.variable_number
            );

            adapter.print_tags();

            if(V.mesh().spatial_dimension() == 2) {
                run_contact<2>(params.contact_params, adapter);
            } else if(V.mesh().spatial_dimension() == 3) {
                run_contact<3>(params.contact_params, adapter);
            }

        });
        
    }
}


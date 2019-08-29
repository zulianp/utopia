#include "utopia_ElasticityApp.hpp"


#include "utopia_ContactSolver.hpp"


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
#include "utopia_ContactStress.hpp"

#include "libmesh/mesh_refinement.h"

namespace utopia {

    class ElasticityAppInput : public Configurable {
    public:
        using ProductSpaceT    = utopia::ProductFunctionSpace<LibMeshFunctionSpace>;
        using MaterialT        = utopia::UIMaterial<ProductSpaceT, USparseMatrix, UVector>;
        using ForcingFunctionT = utopia::UIForcingFunction<ProductSpaceT, UVector>;

        ElasticityAppInput(libMesh::Parallel::Communicator &comm) :
        max_it(20),
        mesh_(comm),
        space_(make_ref(mesh_))
        {}

        void read(Input &is) override
        {
            try {

                is.get("mesh", mesh_);
                is.get("space", space_);

                auto model            = make_unique<MaterialT>(space_.space());
                auto forcing_function = make_unique<ForcingFunctionT>(space_.space());

                is.get("model", *model);
                is.get("forcing-functions", *forcing_function);

                assert(model->good());

                model_ = std::make_shared<ForcedMaterial<USparseMatrix, UVector>>(
                    std::move(model),
                    std::move(forcing_function)
                );

                is.get("max-it", max_it);

            } catch(const std::exception &ex) {
                std::cerr << ex.what() << std::endl;
                assert(false);
            }
        }

        inline bool empty() const
        {
            return mesh_.empty();
        }

        inline libMesh::MeshBase &mesh()
        {
            return mesh_.mesh();
        }

        inline ProductSpaceT &space()
        {
            return space_.space();
        }

        std::shared_ptr< ElasticMaterial<USparseMatrix, UVector> > model()
        {
            return model_;
        }

        void describe(std::ostream &os = std::cout) const
        {

        }

        int max_it;
    private:
        UIMesh<libMesh::DistributedMesh> mesh_;
        UIFunctionSpace<LibMeshFunctionSpace> space_;
        std::shared_ptr< ElasticMaterial<USparseMatrix, UVector> > model_;
    };

    void ElasticityApp::run(Input &in)
    {
        ElasticityAppInput sim_in(comm());
        in.get("elasticity", sim_in);
        sim_in.describe(std::cout);

        if(sim_in.empty()) {
            utopia_error("incomplete input");
            return;
        }

        auto &V = sim_in.space();
        auto &Vx = V[0];
        auto &dof_map = Vx.dof_map();
        auto &material = *sim_in.model();

        UVector x, g, c;
        USparseMatrix H;

        x = ghosted(dof_map.n_dofs(), dof_map.n_local_dofs(), dof_map.get_send_list());

        material.assemble_hessian_and_gradient(x, H, g);
        g *= -1.0;
        apply_boundary_conditions(dof_map, H, g);

        Factorization<USparseMatrix, UVector> solver;

        if(!solver.solve(H, g, c)) {
            write("H.m", H);
            write("g.m", g);
            utopia_error("LinearSolver(0) failed!");
            return;
        }

        x += c;

        if(!material.is_linear()) {
            
            for(int i = 0; i < sim_in.max_it; ++i) {

                material.assemble_hessian_and_gradient(x, H, g);
                g *= -1.0;
                
                apply_boundary_conditions(dof_map, H, x);
                apply_zero_boundary_conditions(dof_map, g);

                c.set(0.0);

                if(!solver.solve(H, g, c)) {
                    write("H.m", H);
                    write("g.m", g);
                    std::cerr << "LinearSolver(" << (i + 1) << ") failed!" << std::endl;;
                    return;
                }

                x += c;

                double norm_g = norm2(g);

                std::cout << "norm_g(" << i << "): " << norm_g << std::endl;

                if(norm_g < 1e-6) {
                    break;
                }

            }
        }

        write("elasticity.e", Vx, x);
    }
}


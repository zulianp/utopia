#include "utopia_ContactApp.hpp"


#include "utopia_ContactSolver.hpp"

#ifndef WITH_TRILINOS_ALGEBRA

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




    class SimulationInput : public Configurable {
    public:
        using ProductSpaceT    = utopia::ProductFunctionSpace<LibMeshFunctionSpace>;
        using MaterialT        = utopia::UIMaterial<ProductSpaceT, USparseMatrix, UVector>;
        using ForcingFunctionT = utopia::UIForcingFunction<ProductSpaceT, UVector>;
        using ContactStressT   = utopia::ContactStress<ProductFunctionSpace<LibMeshFunctionSpace>, USparseMatrix, UVector>;

        SimulationInput(libMesh::Parallel::Communicator &comm) :
        mesh_(comm),
        space_(make_ref(mesh_)),
        dt_(0.1),
        use_amg_(false),
        use_newton(false),
        export_results(true),
        is_steady(false),
        use_gmres(false) {}

        void read(Input &is) override
        {
            try {

                is.get("mesh", mesh_);
                is.get("space", space_);
                is.get("contact", params_);

                auto model            = make_unique<MaterialT>(space_.space());
                auto forcing_function = make_unique<ForcingFunctionT>(space_.space());

                is.get("model", *model);
                is.get("forcing-functions", *forcing_function);
                is.get("dt", dt_);
                is.get("use-amg", use_amg_);
                is.get("export", export_results);
                is.get("use-ssnewton", use_newton);
                is.get("is-steady", is_steady);
                is.get("use-gmres", use_gmres);

                is.get("contact-stress", [this, &model](Input &in) {
                    contact_stress_ = std::make_shared<ContactStressT>(space_.space(), model->params);
                });

                assert(model->good());

                model_ = std::make_shared<ForcedMaterial<USparseMatrix, UVector>>(
                    std::move(model),
                    std::move(forcing_function)
                );

                is.get("qp-solver", [this](Input &in) {
                    this->qp_solver = std::make_shared<PolymorphicQPSolver<USparseMatrix, UVector>>();
                    this->qp_solver->read(in);
                });

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

        inline const UIContactParams &params() const
        {
            return params_;
        }

        std::shared_ptr< ElasticMaterial<USparseMatrix, UVector> > model()
        {
            return model_;
        }

        double dt() const
        {
            return dt_;
        }

        inline bool use_amg() const
        {
            return use_amg_;
        }

        void describe(std::ostream &os = std::cout) const
        {
            // mesh_.describe(os);
            // space_.describe(os);
            os << dt_ << std::endl;
            os << use_amg_ << std::endl;
            os << use_newton << std::endl;
            os << export_results << std::endl;
            os << is_steady << std::endl;
        }

        std::shared_ptr< ContactStressT > contact_stress_;
    private:
        UIMesh<libMesh::DistributedMesh> mesh_;
        UIFunctionSpace<LibMeshFunctionSpace> space_;
        UIContactParams params_;
        std::shared_ptr< ElasticMaterial<USparseMatrix, UVector> > model_;
        
        double dt_;
        bool use_amg_;
    public:
        bool use_newton;
        bool export_results;
        bool is_steady;
        bool use_gmres;

        std::shared_ptr< QPSolver<USparseMatrix, UVector> > qp_solver;
    };


    static void solve_steady(SimulationInput &sim_in)
    {
        typedef utopia::ContactSolver<USparseMatrix, UVector> ContactSolverT;

        const auto &params = sim_in.params();

        auto qp_solver = sim_in.qp_solver;

        ContactSolverT sc(
            make_ref(sim_in.space()),
            sim_in.model(),
            params.contact_params
        );

        if(sim_in.contact_stress_) {
            sc.set_contact_stress(sim_in.contact_stress_);
        }

        if(sim_in.export_results) {
            sc.export_results(true);
        }

        if(qp_solver) {
            sc.set_qp_solver(qp_solver);
        }

        sc.set_tol(1e-3);
        sc.set_max_outer_loops(30);
        sc.set_use_ssn(sim_in.use_newton);

        if(sim_in.use_gmres) {
            // sc.tao().set_ksp_types("gmres", "bjacobi", "petsc");
            sc.set_linear_solver(std::make_shared<GMRES<USparseMatrix, UVector>>("bjacobi"));
        }

#ifdef WITH_M3ELINSOL

        if(sim_in.use_amg()) {
            auto ls = std::make_shared<ASPAMG<USparseMatrix, UVector>>();
            ls->verbose(true);
            auto amg_in_ptr = open_istream("../data/contact/amg_settings.xml");

            if(amg_in_ptr) {
                std::cout << "Using settings" << std::endl;
                amg_in_ptr->get("amg", *ls);
            }

            sc.set_linear_solver(ls);
            sc.set_use_ssn(true);
        }

#endif //WITH_M3ELINSOL

        sc.solve_steady();

    }

    static void solve_transient(SimulationInput &sim_in)
    {
        typedef utopia::ContactStabilizedNewmark<USparseMatrix, UVector> ContactSolverT;

        const auto &params = sim_in.params();

        ContactSolverT sc(
            make_ref(sim_in.space()),
            sim_in.model(),
            sim_in.dt(),
            params.contact_params
        );

        if(sim_in.contact_stress_) {
            sc.set_contact_stress(sim_in.contact_stress_);
        }

        if(sim_in.export_results) {
            sc.export_results(true);
        }

        sc.set_tol(1e-3);
        sc.set_max_outer_loops(30);
        sc.set_use_ssn(sim_in.use_newton);

        if(sim_in.use_gmres) {
            // sc.tao().set_ksp_types("gmres", "bjacobi", "petsc");
            sc.set_linear_solver(std::make_shared<GMRES<USparseMatrix, UVector>>("bjacobi"));
        }

#ifdef WITH_M3ELINSOL

        if(sim_in.use_amg()) {
            auto ls = std::make_shared<ASPAMG<USparseMatrix, UVector>>();
            ls->verbose(true);
            auto amg_in_ptr = open_istream("../data/contact/amg_settings.xml");

            if(amg_in_ptr) {
                std::cout << "Using settings" << std::endl;
                amg_in_ptr->get("amg", *ls);
            }

            sc.set_linear_solver(ls);
            sc.set_use_ssn(true);
        }

#endif //WITH_M3ELINSOL

        sc.initial_condition(1.);
        sc.solve_dynamic(params.n_transient_steps);
    }

    void ContactApp::run(Input &in)
    {
        SimulationInput sim_in(comm());
        in.get("contact-problem", sim_in);
        sim_in.describe(std::cout);

        if(sim_in.is_steady) {
            solve_steady(sim_in);
        } else {
            solve_transient(sim_in);
        }
    }
}

#else

namespace utopia {

    void ContactApp::run(const std::string &path)
    {
        std::cerr << "DOING nothing for trilinos algebra" << std::endl;
    }
}


#endif //WITH_TRILINOS_ALGEBRA


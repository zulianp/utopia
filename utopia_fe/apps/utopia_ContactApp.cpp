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

#include "libmesh/mesh_refinement.h"

namespace utopia {

	typedef utopia::ContactStabilizedNewmark<USparseMatrix, UVector> ContactSolverT;


	class SimulationInput : public Configurable {
	public:
		using ProductSpaceT    = utopia::ProductFunctionSpace<LibMeshFunctionSpace>;
		using MaterialT        = utopia::UIMaterial<ProductSpaceT, USparseMatrix, UVector>;
		using ForcingFunctionT = UIForcingFunction<ProductSpaceT, UVector>;

		SimulationInput(libMesh::Parallel::Communicator &comm) : mesh_(comm), space_(make_ref(mesh_)), dt_(0.1), use_amg_(false), use_newton(false), export_results(true) {}

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

		        model_ = std::make_shared<ForcedMaterial<USparseMatrix, UVector>>(
		            std::move(model),
		            std::move(forcing_function)
		        );

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
	};

	void ContactApp::init(libMesh::LibMeshInit &init)
	{
	    comm_ = make_ref(init.comm());
	}

	void ContactApp::run(const std::string &path)
	{
		SimulationInput sim_in(*comm_);
		auto in_ptr = open_istream(path);//"../data/contact/default.xml");
		in_ptr->get("contact-problem", sim_in);

		const auto &params = sim_in.params();
		
		ContactSolverT sc(
			make_ref(sim_in.space()),
			sim_in.model(),
			sim_in.dt(),
			params.contact_params
		);

		if(sim_in.export_results) {
			sc.export_results(true);
		}

		sc.set_tol(1e-3);
		sc.set_max_outer_loops(30);
		sc.set_use_ssn(sim_in.use_newton);

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
}

#else

namespace utopia {

	void ContactApp::run(const std::string &path)
	{
		std::cerr << "DOING nothing for trilinos algebra" << std::endl;
	}
}


#endif //WITH_TRILINOS_ALGEBRA


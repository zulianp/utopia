#include "utopia_FractureFlow.hpp"

namespace utopia {

	FractureFlow::FractureFlow(libMesh::Parallel::Communicator &comm)
	: mesh(comm), space(make_ref(mesh))
	{}

	void FractureFlow::read(Input &is)
	{
		try {
			is.get("mesh", mesh);
			is.get("space", space);

			auto grid_sampler = std::make_shared<UIScalarSampler<double>>();
			is.get("sampler", *grid_sampler);


			if(!grid_sampler->empty()) {
				sampler = grid_sampler;
			} else {
				auto subdomain_fun = utopia::make_unique<UISubdomainFunction<double>>();

				is.get("diffusivity-blocks", *subdomain_fun);

				if(!subdomain_fun->good()) {
					sampler = std::make_shared<UIConstantFunction<double>>(1.);
				} else {

					if(!subdomain_fun->has_default()) {
						subdomain_fun->set_default(utopia::make_unique<UIConstantFunction<double>>(1.));
					}

					sampler = std::move(subdomain_fun);
				}
			}

			forcing_function = std::make_shared< UIForcingFunction<FunctionSpaceT, UVector> >(space.subspace(0));
			is.get("forcing-function", *forcing_function);

	            //material parameters
			double diffusivity = 1.;
			double diffusivities[3] = {1., 1., 1.};

			is.get("diffusivity", diffusivity);
			is.get("diffusivity-x", diffusivities[0]);
			is.get("diffusivity-y", diffusivities[1]);
			is.get("diffusivity-z", diffusivities[2]);

			int dim = space.subspace(0).mesh().spatial_dimension();

			diffusion_tensor = identity(dim, dim);

			{
				Write<ElementMatrix> w(diffusion_tensor);
				for(int i = 0; i < dim; ++i) {
					diffusion_tensor.set(i, i, diffusivities[i] * diffusivity);
				}
			}

		} catch(const std::exception &ex) {
			std::cerr << ex.what() << std::endl;
			assert(false);
		}
	}

	bool FractureFlow::empty() const
	{
		return mesh.empty();
	}

	void FractureFlow::describe(std::ostream &os) const
	{
		os << "-----------------------------------\n";
	        // mesh.describe(os);
	        // space.describe(os);
	        // forcing_function.describe(os);
		os << "permeability: " << std::endl;

		{
			int dim = size(diffusion_tensor).get(0);
			Read<ElementMatrix> w(diffusion_tensor);
			for(int i = 0; i < dim; ++i) {
				os << diffusion_tensor.get(i, i) << " ";
			}
		}

		os << "\n";
		os << "-----------------------------------\n";
	}

	FractureFlowAuxSystem::FractureFlowAuxSystem(LibMeshFunctionSpace &V, const std::string &name)
	: aux_( V.equation_systems().add_system<libMesh::LinearImplicitSystem>("aux") )
	{
	    var_nums_.push_back( aux_.add_variable(name, libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
	    aux_.init();
	}

	void FractureFlowAuxSystem::sample(const std::shared_ptr<UIFunction<double>> &sampler)
	{
	    LibMeshFunctionSpace V_aperture(aux_, var_nums_[0]);
	    V_aperture.initialize();

	    auto &dof_map = V_aperture.dof_map();

	    sampled = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
	    // sampled = local_zeros(dof_map.n_local_dofs());

	    auto constant_sampler = std::dynamic_pointer_cast<UIConstantFunction<double>>(sampler);

	    if(constant_sampler) {
	        sampled.set(constant_sampler->value());
	    } else {
	        auto u = trial(V_aperture);
	        auto v = test(V_aperture);

	        auto lform = inner(ctx_fun(sampler), v) * dX;

	        UVector aperture_h = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
	        utopia::assemble(lform, aperture_h);

	        USparseMatrix mass_mat;
	        utopia::assemble(inner(u, v) * dX, mass_mat);
	        UVector d_inv = 1./sum(mass_mat, 1);
	        sampled = e_mul(d_inv, aperture_h);
	    }

	    utopia::convert(sampled, *aux_.solution);
	    aux_.solution->close();
	}


}

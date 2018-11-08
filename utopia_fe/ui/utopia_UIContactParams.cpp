#include "utopia_UIContactParams.hpp"

namespace utopia {

	void UIContactParams::read(Input &is)
	{
		std::set<int> temp;

		is.get("radius", contact_params.search_radius);

		std::string type;
		is.get("type", type);

		step_tol = 5e-6;
		is.get("step-tol", step_tol);

		max_nl_iter = 30;
		is.get("max-nl-iter", max_nl_iter);

		is_steady = false;
		n_transient_steps = 1;

		if(type == "steady") {
			is_steady = true;
		}

		use_pg = false;
		std::string solver;
		is.get("solver", solver);

		if(solver == "pg") {
			use_pg = true;
		}

		is.get("n-transient-steps", n_transient_steps);

		is.get("pairs", [this,&temp](Input &is) {
			is.get_all([this,&temp](Input &is) {
				int master = -1, slave = -1;
				is.get("master", master);
				is.get("slave", slave);

		                    // std::cout << master << " " << slave << std::endl;

				assert(master != -1);
				assert(slave  != -1);
				temp.insert(master);
				temp.insert(slave);

				contact_params.contact_pair_tags.push_back({ master, slave });
			});
		});


		contact_surfaces.clear();
		contact_surfaces.insert(contact_surfaces.end(), temp.begin(), temp.end());
	}
}
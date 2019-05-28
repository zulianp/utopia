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

        contact_params.is_glue = std::make_shared<moonolith::IsGlue>();

        is.get("n-transient-steps", n_transient_steps);

        is.get("pairs", [this,&temp](Input &is) {
            is.get_all([this,&temp](Input &is) {
                int master = -1, slave = -1;
                is.get("master", master);
                is.get("slave", slave);

                bool is_glued = false;
                is.get("glue", is_glued);

                assert(master != -1);
                assert(slave  != -1);
                temp.insert(master);
                temp.insert(slave);

                contact_params.contact_pair_tags.push_back({ master, slave });
                contact_params.glued.push_back(is_glued);
                contact_params.is_glue->insert(master, slave);
            });
        });

        contact_surfaces.clear();
        contact_surfaces.insert(contact_surfaces.end(), temp.begin(), temp.end());

        is.get("search-radius", [this](Input &in) {
            in.get("default", contact_params.search_radius);

            contact_params.side_set_search_radius = std::make_shared<moonolith::SearchRadius<double>>(contact_params.search_radius);

            in.get("sides", [this](Input &is) {
                is.get_all([this](Input &is) {
                    int id = -1;
                    double value = contact_params.search_radius;

                    is.get("id", id);
                    is.get("value", value);

                    contact_params.side_set_search_radius->insert(id, value);
                });
            });
        });

        if(!contact_params.side_set_search_radius) {
            contact_params.side_set_search_radius = std::make_shared<moonolith::SearchRadius<double>>(contact_params.search_radius);
        }
    }
}

#ifndef UTOPIA_UI_CONTACT_PARAMS_HPP
#define UTOPIA_UI_CONTACT_PARAMS_HPP

#include "utopia_ui.hpp"
#include "utopia_Contact.hpp"

#include <vector>
#include <iostream>

namespace utopia {

    class UIContactParams final : public Configurable {
    public:
        void read(Input &is) override;
        void describe(std::ostream &os = std::cout) const;

        UIContactParams();

        ContactParams contact_params;
        std::vector<int> contact_surfaces;
        bool is_steady;
        int n_transient_steps;
        double step_tol;
        int max_nl_iter;
        bool use_pg;
    };

}

#endif //UTOPIA_UI_CONTACT_PARAMS_HPP


#include "utopia_IContact.hpp"

namespace utopia {

    void ContactParams::describe(std::ostream &os) const {
        os << "search_radius: " << search_radius << "\n";
        os << "variable_number: " << variable_number << "\n";
        os << "use_biorthogonal_basis: " << use_biorthogonal_basis << "\n";
        os << "master, slave:\n";
        for (const auto &p : contact_pair_tags) {
            os << p.first << ", " << p.second << "\n";
        }

        if (side_set_search_radius) {
            os << "search-radius: " << std::endl;
            side_set_search_radius->describe(os);
        }

        os << std::endl;
    }

}  // namespace utopia

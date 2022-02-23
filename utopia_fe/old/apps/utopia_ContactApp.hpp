#ifndef UTOPIA_CONTACT_APP
#define UTOPIA_CONTACT_APP

#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {

    class ContactApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "-contact"; }
    };
}  // namespace utopia

#endif  // UTOPIA_CONTACT_APP

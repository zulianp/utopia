#ifndef UTOPIA_SMOOTH_APP
#define UTOPIA_SMOOTH_APP

#include "utopia_FEApp.hpp"
#include <string>


namespace utopia {

    class SmoothApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "-smooth";
        }
    };
}


#endif //UTOPIA_SMOOTH_APP

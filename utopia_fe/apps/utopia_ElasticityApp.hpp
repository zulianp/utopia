#ifndef UTOPIA_ELASTICITY_APP
#define UTOPIA_ELASTICITY_APP

#include "utopia_FEApp.hpp"
#include <string>


namespace utopia {

    class ElasticityApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "-elasticity";
        }
    };
}


#endif //UTOPIA_ELASTICITY_APP

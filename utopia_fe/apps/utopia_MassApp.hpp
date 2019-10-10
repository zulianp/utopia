#ifndef UTOPIA_MASS_APP_HPP
#define UTOPIA_MASS_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"
#include "utopia_UIForcingFunction.hpp"


namespace utopia {
    class MassApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "-mass";
        }   
    };




    

   

}


#endif //UTOPIA_MASS_APP_HPP

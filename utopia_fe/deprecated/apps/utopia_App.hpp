#ifndef UTOPIA_APP_HPP
#define UTOPIA_APP_HPP

#include <string>
#include "utopia_ui.hpp"

namespace utopia {
    class App {
    public:
        virtual ~App() {}
        virtual void run(Input &in) = 0;
    };
}  // namespace utopia

#endif  // UTOPIA_APP_HPP

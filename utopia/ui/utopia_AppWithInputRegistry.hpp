#ifndef UTOPIA_APP_WITH_INPUT_REGISTRY_H
#define UTOPIA_APP_WITH_INPUT_REGISTRY_H

#include <map>
#include <string>
#include <functional>
#include <iostream>

#include "utopia_ActionRegistry.hpp"
#include "utopia_NaryActionRegistry.hpp"
#include "utopia_Input.hpp"

namespace utopia {
    using AppWithInputRegistry = utopia::NaryActionRegistry<Input &>;
}

#endif //UTOPIA_APP_WITH_INPUT_REGISTRY_H


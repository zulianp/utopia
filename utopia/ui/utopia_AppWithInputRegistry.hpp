#ifndef UTOPIA_APP_WITH_INPUT_REGISTRY_H
#define UTOPIA_APP_WITH_INPUT_REGISTRY_H

#include <functional>
#include <iostream>
#include <map>
#include <string>

#include "utopia_ActionRegistry.hpp"
#include "utopia_Input.hpp"
#include "utopia_NaryActionRegistry.hpp"

namespace utopia {
    using AppWithInputRegistry = utopia::NaryActionRegistry<Input &>;
}

#endif  // UTOPIA_APP_WITH_INPUT_REGISTRY_H

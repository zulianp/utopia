#include "../utopia.hpp"
#include "utopia_Instance.hpp"
#include "utopia_Layout.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

#include <iostream>

namespace scripting {
    Layout::Layout(class Communicator, const int n_global, const int n_local) : impl_(nullptr) {}

    Layout::~Layout() { delete impl_; }

}  // namespace scripting
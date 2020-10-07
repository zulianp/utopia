#include "../utopia.hpp"
#include "utopia_Instance.hpp"
#include "utopia_Layout.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

#include <iostream>

namespace scripting {
    Layout::Layout() : impl_(nullptr) {}

    Layout::~Layout() { delete impl_; }

    // SerialLayout SerialLayout::serial_layout(const SizeType &size) { return impl_->serial_layout(size); }

}  // namespace scripting
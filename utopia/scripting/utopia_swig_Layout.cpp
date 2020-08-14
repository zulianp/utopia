#include "../utopia.hpp"
#include "utopia_Instance.hpp"
#include "utopia_Layout.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

#include <iostream>

namespace scripting {
    Layout::Layout() : impl_(nullptr) {
        auto layout = Factory::new_layout();
        if (!layout) {
            utopia::out() << "[Error] Layout could not be constructed" << std::endl;
            return;
        }

        impl_ = layout.get();

        layout.release();
    }

    Layout::~Layout() { delete impl_; }

    // Layout::serial_layout(const SizeType &size) { return impl_->serial_layout(size); }
}  // namespace scripting
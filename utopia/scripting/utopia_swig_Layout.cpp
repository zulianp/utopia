#include "../utopia.hpp"
#include "utopia_Instance.hpp"
#include "utopia_Layout.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

#include <iostream>

namespace scripting {
    Layout::Layout() : impl_(nullptr) {
        auto comm = Factory::new_layout();
        if (!comm) {
            utopia::out() << "[Error] Communicator could not be constructed" << std::endl;
            return;
        }

        impl_ = comm.get();

        comm.release();
    }

    Layout::~Layout() { delete impl_; }

    // Layout::serial_layout(const SizeType &size) { return impl_->serial_layout(size); }
}  // namespace scripting
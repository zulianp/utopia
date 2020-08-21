#include "../utopia.hpp"
#include "utopia_Communicator.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_script.hpp"

#include <iostream>

namespace scripting {
    Communicator::Communicator() : impl_(nullptr) {
        auto comm = Factory::new_communicator();
        if (!comm) {
            utopia::out() << "[Error] Communicator could not be constructed" << std::endl;
            return;
        }

        impl_ = comm.get();

        comm.release();
    }

    Communicator::~Communicator() { delete impl_; }

    int Communicator::rank() { return impl_->rank(); }
    int Communicator::size() { return impl_->size(); }

    SelfCommunicator::SelfCommunicator() : impl_(nullptr) {
        auto selfComm = Factory::new_communicator();
        if (!selfComm) {
            utopia::out() << "[Error] SelfCommunicator could not be constructed" << std::endl;
            return;
        }

        // impl_ = selfComm.get();

        selfComm.release();
    }

    SelfCommunicator::~SelfCommunicator() { delete impl_; }

    // Communicator SelfCommunicator::get_default() { return impl_->get_default(); }

    int SelfCommunicator::rank() { return impl_->rank(); }
    int SelfCommunicator::size() { return impl_->size(); }

    // void Communicator::get_the_default() { return impl_->get_default(); }
}  // namespace scripting
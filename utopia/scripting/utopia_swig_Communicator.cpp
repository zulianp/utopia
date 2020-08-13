#include "../utopia.hpp"
#include "utopia_Communicator.hpp"
#include "utopia_Instance.hpp"
#include "utopia_Layout.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

#include <iostream>

namespace utopia {

#ifdef WITH_PETSC
    UTOPIA_FACTORY_REGISTER_VECTOR(PetscVector);
    UTOPIA_FACTORY_REGISTER_MATRIX(PetscMatrix);
#endif  // WITH_PETSC

#ifdef WITH_TRILINOS
    UTOPIA_FACTORY_REGISTER_VECTOR(TpetraVector);
#endif  // WITH_PETSC

}  // namespace utopia
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

    void Communicator::get_the_default() { return impl_->get_default(); }
}  // namespace scripting
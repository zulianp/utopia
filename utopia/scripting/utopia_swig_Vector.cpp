#include <iostream>
#include "../utopia.hpp"
#include "utopia_AbstractVector.hpp"
#include "utopia_Instance.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

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
    Vector::Vector() : impl_(nullptr) {
        auto vec = Factory::new_vector();

        if (!vec) {
            utopia::out() << "[Error] Vector could not be constructed" << std::endl;
            return;
        }

        impl_ = vec.get();

        vec.release();
    }

    Vector::~Vector() { delete impl_; }

    void Vector::print_info() { utopia::out() << "Vector::print()" << std::endl; }

    // TODO:: change so that we do not take size as parameter but a layout.
    void Vector::create_vector(const SizeType &size, const Scalar &value) {
        auto l = utopia::serial_layout(size);
        impl_->values(l, value);
    }

    void Vector::clear_vector() { impl_->clear(); }

    bool Vector::is_empty() { return impl_->empty(); }

    void Vector::disp() { impl_->describe(); }

}  // namespace scripting
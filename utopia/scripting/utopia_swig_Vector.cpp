#include <iostream>
#include "../utopia.hpp"
#include "utopia_AbstractVector.hpp"
#include "utopia_Instance.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

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

    // TODO:: change so that we do not take size as parameter but a layout.
    void Vector::create_vector(const SizeType &size, const Scalar &value) {
        auto l = utopia::serial_layout(size);
        impl_->values(l, value);
    }

    void Vector::clear() { impl_->clear(); }

    bool Vector::empty() { return impl_->empty(); }

    void Vector::disp() { impl_->describe(); }

    void Vector::set(const SizeType &i, const Scalar &value) { impl_->set(i, value); }

    Scalar Vector::get(const SizeType &i) { return impl_->get(i); }

    void Vector::describe() { impl_->describe(); }

    // Communicator &Vector::comm() { return impl_->& comm(); }

    Scalar Vector::norm1() { return impl_->norm1(); }

    // void Vector::values(const Layout &l, const Scalar &value) { impl_->values(l, value); }

}  // namespace scripting
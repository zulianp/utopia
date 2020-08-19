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

    SizeType Vector::local_size() { return impl_->local_size(); }

    void Vector::c_set(const SizeType &i, const Scalar &value) { impl_->c_set(i, value); }
    void Vector::c_add(const SizeType &i, const Scalar &value) { impl_->c_add(i, value); }

    void Vector::read_lock() { impl_->read_lock(); }
    void Vector::read_unlock() { impl_->read_unlock(); }

    void Vector::clear() { impl_->clear(); }

    bool Vector::empty() { return impl_->empty(); }

    void Vector::disp() { impl_->describe(); }

    void Vector::set(const SizeType &i, const Scalar &value) { impl_->set(i, value); }

    Scalar Vector::get(const SizeType &i) { return impl_->get(i); }

    void Vector::describe() { impl_->describe(); }

    // Communicator &Vector::comm() { return impl_->& comm(); }

    Scalar Vector::norm_infty() {return impl_->norm_infty();}

    Scalar Vector::norm1() { return impl_->norm1(); }

    // void Vector::e_mul(const Scalar &other) override {impl_->e_mul(other);}

    // void Vector::e_min(const Scalar &other) override {impl_->e_min(other);}

    // void Vector::e_div(const Scalar &other) override {impl_->e_div(other);}

    // void Vector::e_max(const Scalar &other) override {impl_->e_max(other);}

    // void Vector::values(const Layout &l, const Scalar &value) { impl_->values(l, value); }

}  // namespace scripting
#include <iostream>
#include "../utopia.hpp"
#include "utopia_AbstractVector.hpp"
#include "utopia_AbstractMatrix.hpp"
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

    Scalar Vector::norm_infty() { return impl_->norm_infty(); }

    Scalar Vector::norm1() { return impl_->norm1(); }

    // void Vector::e_mul(const Scalar &other) override {impl_->e_mul(other);}

    // void Vector::e_min(const Scalar &other) override {impl_->e_min(other);}

    // void Vector::e_div(const Scalar &other) override {impl_->e_div(other);}

    // void Vector::e_max(const Scalar &other) override {impl_->e_max(other);}

    // void Vector::values(const Layout &l, const Scalar &value) { impl_->values(l, value); }

}  // namespace scripting

namespace scripting {
    Layout::Layout() : impl_(nullptr) {}

    Layout::~Layout() { delete impl_; }

    // SerialLayout SerialLayout::serial_layout(const SizeType &size) { return impl_->serial_layout(size); }

}  // namespace scripting

// namespace utopia
namespace scripting {
    SparseMatrix::SparseMatrix() : impl_(nullptr) {
        auto mat = Factory::new_matrix();

        if (!mat) {
            utopia::out() << "[Error] Matrix could not be constructed" << std::endl;
            return;
        }

        impl_ = mat.get();
        mat.release();
    }

    SparseMatrix::~SparseMatrix() { delete impl_; }

    void SparseMatrix::print_info() { utopia::out() << "SparseMatrix::print()" << std::endl; }

    // // TODO:: Change so that we have a layout as paramater not size.
    // void SparseMatrix::create_identity_matrix(const SizeType &size) {
    //     auto l = utopia::serial_layout(size);
    //     impl_->identity(utopia::square_matrix_layout(l), 1);
    // }

    void SparseMatrix::clear() { impl_->clear(); }

    bool SparseMatrix::empty() { return impl_->empty(); }

    void SparseMatrix::disp() { impl_->describe(); }

    void SparseMatrix::set(const SizeType &i, const SizeType &j, const Scalar &value) { impl_->set(i, j, value); }
    void SparseMatrix::add(const SizeType &i, const SizeType &j, const Scalar &value) { impl_->add(i, j, value); }

    void SparseMatrix::c_set(const SizeType &i, const SizeType &j, const Scalar &value) { impl_->c_set(i, j, value); }

    void SparseMatrix::c_add(const SizeType &i, const SizeType &j, const Scalar &value) { impl_->c_add(i, j, value); }

    SizeType SparseMatrix::cols() { return impl_->cols(); }

    SizeType SparseMatrix::rows() { return impl_->rows(); }

    SizeType SparseMatrix::local_cols() { return impl_->local_cols(); }

    SizeType SparseMatrix::local_rows() { return impl_->local_rows(); }

    void SparseMatrix::read_lock() { impl_->read_lock(); }
    void SparseMatrix::read_unlock() { impl_->read_unlock(); }

    Scalar SparseMatrix::norm_infty() { return impl_->norm_infty(); }

    Scalar SparseMatrix::norm1() { return impl_->norm1(); }
    Scalar SparseMatrix::norm2() { return impl_->norm2(); }

    Scalar SparseMatrix::trace() { return impl_->trace(); }

    // void transform(const Sqrt &op) { return impl_->transform(op); }
}  // namespace scripting
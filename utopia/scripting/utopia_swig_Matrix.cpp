#include <iostream>
#include "../utopia.hpp"
#include "utopia_AbstractMatrix.hpp"
#include "utopia_Instance.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

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
    void SparseMatrix::add(const SizeType &i, const SizeType &j, const Scalar &value){ impl_->add(i, j, value); }

    void SparseMatrix::c_set(const SizeType &i, const SizeType &j, const Scalar &value) { impl_->c_set(i, j, value); }

    void SparseMatrix::c_add(const SizeType &i, const SizeType &j, const Scalar &value){ impl_->c_add(i, j, value); }

    SizeType SparseMatrix::cols() {return impl_->cols();}

    SizeType SparseMatrix::rows() {return impl_->rows();}

    SizeType SparseMatrix::local_cols() {return impl_->local_cols();}

    SizeType SparseMatrix::local_rows() {return impl_->local_rows();}

    void SparseMatrix::read_lock() { impl_->read_lock(); }
    void SparseMatrix::read_unlock() { impl_->read_unlock(); }


    Scalar SparseMatrix::norm_infty() {return impl_->norm_infty();}

    Scalar SparseMatrix::norm1() { return impl_->norm1(); }

}  // namespace scripting
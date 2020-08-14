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

    // TODO:: Change so that we have a layout as paramater not size.
    void SparseMatrix::create_identity_matrix(const SizeType &size) {
        auto l = utopia::serial_layout(size);
        impl_->identity(utopia::square_matrix_layout(l), 1);
    }

    void SparseMatrix::clear_matrix() { impl_->clear(); }

    bool SparseMatrix::is_empty() { return impl_->empty(); }

    void SparseMatrix::disp() { impl_->describe(); }
}  // namespace scripting
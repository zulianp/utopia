#include "utopia_script.hpp"
#include <iostream>
#include "../utopia.hpp"
#include "utopia.hpp"
#include "utopia_AbstractVector.hpp"
#include "utopia_Instance.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

void print_array(double *seq, int n) {
    printf("array with length %d:\n", n);

    for (int i = 0; i < n; ++i) {
        printf("%g\n", seq[i]);
    }
}

namespace utopia {

#ifdef UTOPIA_ENABLE_PETSC
    UTOPIA_FACTORY_REGISTER_VECTOR(PetscVector);
    UTOPIA_FACTORY_REGISTER_MATRIX(PetscMatrix);
#endif  // UTOPIA_ENABLE_PETSC

#ifdef UTOPIA_ENABLE_TRILINOS
    UTOPIA_FACTORY_REGISTER_VECTOR(TpetraVector);
#endif  // UTOPIA_ENABLE_PETSC

}  // namespace utopia

namespace scripting {

    void init(int argc, char *argv[]) {
        using namespace utopia;
        Utopia::Init(argc, argv);
        scripting::print_info();

        scripting::Factory::print_info();
    }

    void init() {
        // FIXME?
        int argc = 1;
        std::string argv = "utopia_script";
        char *argv_ptr = &argv[0];
        init(argc, &argv_ptr);
    }

    void print_info() { utopia::out() << "Utopia\nversion: " << UTOPIA_VERSION << std::endl; }

    void finalize() { utopia::Utopia::Finalize(); }

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

    Layout::Layout(const Communicator &comm, LocalSizeType local_size, SizeType global_size)
        :

          impl_(nullptr) {
        auto layout = std::make_unique<LayoutImpl>(*comm.get_communicator(), local_size, global_size);

        if (!layout) {
            utopia::out() << "[Error] Vector could not be constructed" << std::endl;
            return;
        }

        impl_ = layout.get();
        layout.release();
    }

    Layout::~Layout() { delete impl_; }

    // input inside here, vcector of double
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
    void Vector::values(const Layout &l, const Scalar &value) {
        auto ll = *l.get_layout();
        impl_->values(ll, value);
    }
    void Vector::add(const SizeType &i, const Scalar &value) { impl_->add(i, value); }
    void Vector::axpy(Scalar alpha, Vector *x) { impl_->axpy(alpha, *x->impl_); }
    void Vector::describe() const { impl_->describe(); }
    bool Vector::equals(const Vector *other, const Scalar tol) const { return impl_->equals(*other->impl_, tol); }
    Scalar Vector::dot(const Vector *x) const { return impl_->dot(*x->impl_); }
    void Vector::set(const SizeType &i, const Scalar &value) { impl_->set(i, value); }
    void Vector::convert_into_uvector(Scalar *values, const Layout &l) {
        if (impl_->empty()) {
            impl_->values(*l.get_layout(), 0.0);
        }

        {
            impl_->write_lock(utopia::LOCAL);

            utopia::Range rr = impl_->range();

            for (auto i = rr.begin(); i < rr.end(); ++i) {
                impl_->set(i, *values);
                ++values;
            }

            impl_->write_unlock(utopia::LOCAL);
        }
    }

}  // namespace scripting

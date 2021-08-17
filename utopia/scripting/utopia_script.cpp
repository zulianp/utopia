#include "utopia_script.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "../utopia.hpp"
#include "utopia.hpp"
#include "utopia_AbstractVector.hpp"
#include "utopia_Instance.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

#include <iostream>

void print_array(double *seq, int n) {
    printf("array with length %d :\n", n);
    for (int i = 0; i < n; ++i) {
        printf("%g\n", seq[i]);
    }
}

namespace utopia {

#ifdef UTOPIA_WITH_PETSC
    UTOPIA_FACTORY_REGISTER_VECTOR(PetscVector);
    UTOPIA_FACTORY_REGISTER_MATRIX(PetscMatrix);
#endif  // UTOPIA_WITH_PETSC

#ifdef UTOPIA_WITH_TRILINOS
    UTOPIA_FACTORY_REGISTER_VECTOR(TpetraVector);
#endif  // UTOPIA_WITH_PETSC

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
    Scalar Vector::get(const SizeType &i) const { return impl_->get(i); }
    void Vector::scale(const Scalar &a) { impl_->scale(a); }

    void Vector::print_array(double *seq, int n) {
        printf("array with length %d :\n", n);
        for (int i = 0; i < n; ++i) {
            printf("%g\n", seq[i]);
        }
    }

    SizeType Vector::local_size() { return impl_->local_size(); }

    void Vector::serial_uconversion(double *seq, int n) {
        {
            if (impl_->empty()) {
                auto l = utopia::serial_layout(n);
                impl_->values(l, 0);
            }

            impl_->write_lock(utopia::LOCAL);
            for (auto i = 0; i < n; ++i) {
                impl_->set(i, *seq);
                ++seq;
            }

            impl_->write_unlock(utopia::LOCAL);
        }
    }

    // double *Vector::from_utopia_to_carray() {
    //     int size = impl_->size();
    //     double *temp = new double[size];
    //     {
    //         impl_->write_lock(utopia::LOCAL);
    //         utopia::Range rr = impl_->range();
    //         for (auto i = rr.begin(); i < rr.end(); ++i) {
    //             temp[i] = impl_->get(i);
    //         }

    //         impl_->write_unlock(utopia::LOCAL);
    //     }
    //     for (int j = 0; j < size; j++) {
    //         printf("%g\n", temp[j]);
    //     }
    //     return temp;
    // }

    void Vector::parallel_uconversion(float *values, const Layout &l) {
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

    void Vector::write_into_carray(double *double_array) {
        {
            impl_->write_lock(utopia::LOCAL);
            utopia::Range rr = impl_->range();
            for (auto i = rr.begin(); i < rr.end(); ++i) {
                double_array[i] = impl_->get(i);
            }

            impl_->write_unlock(utopia::LOCAL);
        }
    }

    // void Vector::parallel_uconversion(std::vector<double> values, const Layout &l) {
    //     if (impl_->empty()) {
    //         impl_->values(*l.get_layout(), 0.0);
    //     }

    //     {
    //         impl_->write_lock(utopia::LOCAL);
    //         auto it = values.begin();
    //         utopia::Range rr = impl_->range();
    //         for (auto i = rr.begin(); i < rr.end(); ++i) {
    //             impl_->set(i, *it);
    //             ++it;
    //         }

    //         impl_->write_unlock(utopia::LOCAL);
    //     }
    // }

}  // namespace scripting

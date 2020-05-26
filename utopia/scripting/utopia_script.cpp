#include "utopia_script.hpp"
#include "utopia.hpp"
#include "utopia_Instance.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"

#include <iostream>

namespace scripting {

    void init(int argc, char *argv[]) {
        using namespace utopia;
        Utopia::Init(argc, argv);
        scripting::print_info();

#ifdef WITH_PETSC
        UTOPIA_FACTORY_REGISTER_VECTOR(PetscVector);
        UTOPIA_FACTORY_REGISTER_MATRIX(PetscMatrix);
#endif  // WITH_PETSC

#ifdef WITH_TRILINOS
        UTOPIA_FACTORY_REGISTER_VECTOR(TpetraVector);
#endif  // WITH_PETSC

        scripting::Factory::print_info();
    }

    void init() {
        // FIXME?
        int argc = 1;
        std::string argv = "utopia_script";
        char *argv_ptr = &argv[0];
        init(argc, &argv_ptr);
    }

    void print_info() { std::cout << "Utopia\nversion: " << UTOPIA_VERSION << std::endl; }

    void finalize() { utopia::Utopia::Finalize(); }

    SparseMatrix::SparseMatrix() : impl_(nullptr) {
        auto mat = Factory::new_matrix();

        if (!mat) {
            std::cout << "[Error] Matrix could not be constructed" << std::endl;
            return;
        }

        impl_ = mat.get();
        mat.release();
    }

    SparseMatrix::~SparseMatrix() { delete impl_; }

    void SparseMatrix::print_info() { std::cout << "SparseMatrix::print()" << std::endl; }

    Vector::Vector() : impl_(nullptr) {
        auto vec = Factory::new_vector();

        if (!vec) {
            std::cout << "[Error] Vector could not be constructed" << std::endl;
            return;
        }

        impl_ = vec.get();
        vec.release();
    }

    Vector::~Vector() { delete impl_; }

    void Vector::print_info() { std::cout << "Vector::print()" << std::endl; }

}  // namespace scripting

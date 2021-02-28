#include "utopia_script.hpp"
#include <iostream>
#include "../utopia.hpp"
#include "utopia.hpp"
#include "utopia_AbstractVector.hpp"
#include "utopia_Instance.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"
#include "utopia_script.hpp"

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

        if(!comm) {
            utopia::out() << "[Error] Communicator could not be constructed" << std::endl;
            return; 
        }
        impl_ = comm.get();
        comm.release();
    }

    Communicator::~Communicator() {delete impl_;}

    Layout::Layout(const Communicator &comm, int Order, LocalSizeType local_size, SizeType global_size) : 
    
        impl_(nullptr), comm_(comm), Order_(Order), local_size_(local_size), global_size_(global_size) {

            auto layout = std::make_unique<LayoutImpl>(*comm_.get_communicator(), local_size_, global_size_); 

        if (!layout) {
            utopia::out() << "[Error] Vector could not be constructed" << std::endl;
            return;
        }

        impl_ = layout.get();
        layout.release();
    }

    Layout::~Layout() {delete impl_;}

    






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
    void Vector::set(const Scalar &val) { impl_->set(val); }

}  // namespace scripting

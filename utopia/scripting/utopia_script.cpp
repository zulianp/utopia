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

    Layout::Layout(const Communicator *comm_var, int Order_var, LocalSizeType local_size_var, SizeType global_size_var) {
       
    // impl_->utopia::Communicator = *comm_var;
    // impl_->Order = Order_var;
    // impl_->local_size = local_size_var;
    // impl_->global_size = global_size_var;
    //impl_ = utopia::layout(&comm_var, Order_var, local_size_var, global_size_var);




        // auto layout = utopia::layout(*comm_var, Order_var, local_size_var, global_size_var);

        // impl_ = layout.get();
        // layout.release();

        *comm = *comm_var;
        local_size = local_size_var;
        global_size = global_size_var;
        int Order = Order_var;
    }

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

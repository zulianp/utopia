
#ifndef UTOPIA_UTOPIA_PETSCVECTOR_H
#define UTOPIA_UTOPIA_PETSCVECTOR_H

#include "utopia_petsc_Error.hpp"

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"

#include "petscmat.h"
#include "petscvec.h"

#include <memory>

namespace utopia {

    class PETScVector {
    public:
        inline PETScVector()
        : vec_(nullptr), initialized_(false)
        {
        }

        inline ~PETScVector()
        {
            destroy();
        }

        PETScVector(const PETScVector &other)
        {
            if(other.vec_) {
                PETScError::Check(VecDuplicate(other.vec_, &vec_));
                PETScError::Check(VecCopy(other.vec_, vec_));
                initialized_ = other.initialized_;
            } else {
                vec_ = nullptr;
                initialized_ = false;
            }
        }

        inline VecType type() const
        {
            assert(vec_ != nullptr);
            VecType ret;
            VecGetType(vec_, &ret);
            return ret;
        }

        inline PetscInt local_size() const
        {
            PetscInt ret;
            VecGetLocalSize(vec_, &ret);
            return ret;
        }

        inline PetscInt size() const
        {
            PetscInt ret;
            VecGetSize(vec_, &ret);
            return ret;
        }

        inline MPI_Comm communicator() const {
            assert(vec_ != nullptr);
            MPI_Comm comm = PetscObjectComm((PetscObject) implementation());
            assert(comm != MPI_COMM_NULL);
            return comm;
        }

        // assign operator
       inline PETScVector &operator=(const PETScVector &other) {
            if(this == &other) return *this;
            destroy();

            if(other.vec_) {
                PETScError::Check(VecDuplicate(other.vec_, &vec_));
                PETScError::Check(VecCopy(other.vec_, vec_));
            }

            return *this;
        }

        inline PETScVector &operator=(PETScVector &&other) {
            if(this == &other) return *this;
            destroy();

            vec_ = other.vec_;
            other.vec_ = nullptr;
            return *this;
        }

        inline void destroy() {
            if(vec_) {
                VecDestroy(&vec_);
                vec_ = nullptr;
            }

            initialized_ = false;
        }

        inline Vec &implementation() {
            return vec_;
        }

        inline const Vec &implementation() const {
            assert(vec_ != nullptr);
            return vec_;
        }

        inline void describe() const {
            assert(vec_ != nullptr);
            VecView(vec_, PETSC_VIEWER_STDOUT_(communicator()));
        }

        inline bool is_null() const
        {
            return vec_ == nullptr;
        }

        inline bool initialized() const {
            return initialized_;
        }

        inline void set_initialized(const bool val) 
        {
            initialized_ = val;
        }

    private:
        Vec vec_;
        bool initialized_;
    };

    // class PETScVector {
    // public:
    //     inline PETScVector(MPI_Comm comm = PETSC_COMM_WORLD)
    //     : initialized_(false)
    //     {
    //         PETScError::Check(VecCreate(comm, &vec_));
    //     }

    //     inline ~PETScVector()
    //     {
    //         destroy();
    //     }

    //     PETScVector(const PETScVector &other)
    //     {
    //         PETScError::Check(VecDuplicate(other.vec_, &vec_));
    //         PETScError::Check(VecCopy(other.vec_, vec_));
    //         initialized_ = other.initialized_;
    //     }

    //     inline MPI_Comm communicator() const {
    //         MPI_Comm comm = PetscObjectComm((PetscObject) implementation());
    //         assert(comm != MPI_COMM_NULL);
    //         return comm;
    //     }

    //     // assign operator
    //    inline PETScVector &operator=(const PETScVector &other) {
    //         if(this == &other) return *this;
    //         destroy();
    //         PETScError::Check(VecDuplicate(other.vec_, &vec_));
    //         PETScError::Check(VecCopy(other.vec_, vec_));
    //         return *this;
    //     }

    //     inline PETScVector &operator=(PETScVector &&other) {
    //         if(this == &other) return *this;
    //         destroy();
    //         vec_ = other.vec_;
    //         other.vec_ = nullptr;
    //         return *this;
    //     }

    //     inline void destroy() {
    //         VecDestroy(&vec_);
    //         initialized_ = false;
    //     }

    //     inline Vec &implementation() {
    //         return vec_;
    //     }

    //     inline const Vec &implementation() const {
    //         return vec_;
    //     }

    //     inline void describe() const {
    //         VecView(vec_, PETSC_VIEWER_STDOUT_(communicator()));
    //     }

    //     inline bool initialized() const {
    //         return initialized_;
    //     }

    //     inline void set_initialized(const bool val) 
    //     {
    //         initialized_ = val;
    //     }

    // private:
    //     Vec vec_;
    //     bool initialized_;
    // };
}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H

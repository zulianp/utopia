
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
        inline PETScVector(MPI_Comm comm = PETSC_COMM_WORLD)
        : initialized_(false)
        {
            PETScError::Check(VecCreate(comm, &_vec));
        }

        inline ~PETScVector()
        {
            destroy();
        }

        PETScVector(const PETScVector &other)
        {
            PETScError::Check(VecDuplicate(other._vec, &_vec));
            PETScError::Check(VecCopy(other._vec, _vec));
            initialized_ = other.initialized_;
        }

        inline MPI_Comm communicator() const {
            MPI_Comm comm = PetscObjectComm((PetscObject) implementation());
            assert(comm != MPI_COMM_NULL);
            return comm;
        }

        // assign operator
       inline PETScVector &operator=(const PETScVector &other) {
            if(this == &other) return *this;
            destroy();
            PETScError::Check(VecDuplicate(other._vec, &_vec));
            PETScError::Check(VecCopy(other._vec, _vec));
            return *this;
        }

        inline PETScVector &operator=(PETScVector &&other) {
            if(this == &other) return *this;
            destroy();
            _vec = other._vec;
            other._vec = nullptr;
            return *this;
        }

        inline void destroy() {
            VecDestroy(&_vec);
            initialized_ = false;
        }

        inline Vec &implementation() {
            return _vec;
        }

        inline const Vec &implementation() const {
            return _vec;
        }

        inline void describe() const {
            VecView(_vec, PETSC_VIEWER_STDOUT_WORLD);
        }

        inline bool initialized() const {
            return initialized_;
        }

        inline void set_initialized(const bool val) 
        {
            initialized_ = val;
        }

    private:
        Vec _vec;
        bool initialized_;
    };
}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H

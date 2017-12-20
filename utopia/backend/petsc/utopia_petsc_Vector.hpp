
#ifndef UTOPIA_UTOPIA_PETSCVECTOR_H
#define UTOPIA_UTOPIA_PETSCVECTOR_H

#include "utopia_petsc_Error.hpp"

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"

#include "petscmat.h"
#include "petscvec.h"

#include <memory>

namespace utopia {

    // New PETScVector (one level wrapper)

    class PETScVector {

    public:
        PETScVector(MPI_Comm comm = PETSC_COMM_WORLD)
            : _comm(comm), _initialized(false), _nLocalGhosts(0) {
            // FIXME Check error on create????
            PETScError::Check(VecCreate(_comm, &_vec));
        }

        ~PETScVector() {
            destroy();
        }

        PETScVector(const PETScVector &other) {
            // Copy communicator
            _comm = other._comm;
            // Destory current vector
            //PETScError::Check(VecDestroy(&_vec));
            // Copy other
            PETScError::Check(VecDuplicate(other._vec, &_vec));
            PETScError::Check(VecCopy(other._vec, _vec));
            // Copy fields
            // FIXME: Need to copy this fields??
            _initialized = other._initialized;
            _nLocalGhosts = other._nLocalGhosts;
            _ghosts = other._ghosts;
        }

        // get communicator
        MPI_Comm &communicator() {
            return _comm;
        }

        const MPI_Comm &communicator() const {
            return _comm;
        }

        // set communicator for MPI
        void setCommunicator(const MPI_Comm comm) {
            _comm = comm;
        }

        // assign operator
        PETScVector &operator=(const PETScVector &other) {
            if(this == &other) return *this;
            // TODO
            //_wrapper = std::make_shared<PETScVectorWrapper>();
            //other._wrapper->copy(*_wrapper);
            destroy();
            PETScError::Check(VecDuplicate(other._vec, &_vec));
            PETScError::Check(VecCopy(other._vec, _vec));
            _comm = other._comm;

            return *this;
        }

        PETScVector &operator=(PETScVector &&other) {
            if(this == &other) return *this;
            // TODO
            //_wrapper = std::make_shared<PETScVectorWrapper>();
            //other._wrapper->copy(*_wrapper);
            destroy();
            // PETScError::Check(VecDuplicate(other._vec, &_vec));
            // PETScError::Check(VecCopy(other._vec, _vec));
            _comm = other._comm;
            _vec= other._vec;
            other._vec = nullptr;

            return *this;
        }

        // destroy vector
        void destroy() {
            VecDestroy(&_vec);
        }

        void resize(const PetscInt nLocal, const PetscInt nGlobal = PETSC_DETERMINE) {
            destroy();
            // FIXME Error check?
            PETScError::Check(VecCreateMPI(_comm, nLocal, nGlobal, &_vec));
            _initialized = true;
        }

        // assembly vector
        void assemblyBegin() {
            PETScError::Check(VecAssemblyBegin(_vec));
        }

        void assemblyEnd() {
            PETScError::Check(VecAssemblyEnd(_vec));
        }

        Vec &implementation() {
            return _vec;
        }

        const Vec &implementation() const {
            return _vec;
        }

        void describe() const {
            VecView(_vec, PETSC_VIEWER_STDOUT_WORLD);
        }

        // // Check if vector is initialized
        inline bool isInitialized() const {
            return _initialized;
        }


    // private data
    private:
        MPI_Comm _comm;
        Vec _vec;

        bool _initialized;
        PetscInt _nLocalGhosts;
        std::vector<PetscInt> _ghosts;
    };
}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H

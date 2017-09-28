
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

////        PETScVector(const PetscInt nLocal, MPI_Comm comm = PETSC_COMM_WORLD)
////                : _initialized(true), _nLocalGhosts(0)
////        {
////            _wrapper = std::make_shared<PETScVectorWrapper>(comm);
////
////            long nGlobal = nLocal;
////            MPI_Allreduce(MPI_IN_PLACE, &nGlobal, 1, MPI_LONG, MPI_SUM, _wrapper->communicator());
////            VecCreateMPI(comm, nLocal, nGlobal, &implementation());
////        }

        // initialize vector
        bool initialize(const PetscInt nLocal, const PetscInt nGlobal) {
            // FIXME Error check?
            PETScError::Check(VecCreateMPI(_comm, nLocal, nGlobal, &_vec));
            return true;
        }

        const std::vector<PetscInt> &ghosts() const {
            return _ghosts;
        }

        bool initializeWithGhosts(const PetscInt nLocal, const PetscInt nGlobal, const std::vector<PetscInt> &ghosts) {
            // FIXME Error check?
            PETScError::Check(VecCreate(PETSC_COMM_WORLD, &_vec));
            PETScError::Check(VecSetType(_vec, VECMPI));
            PETScError::Check(VecSetSizes(_vec, nLocal, nGlobal));

            if (!ghosts.empty()) {
                VecMPISetGhost(_vec, ghosts.size(), &ghosts[0]);
                _nLocalGhosts = ghosts.size();
                _ghosts = ghosts;
                updateGhosts();
            }

            return true;
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

        void finalize() {
            assemblyBegin();
            assemblyEnd();
        }

        void setGlobalValue(const PetscInt i, const PetscReal value) {
            PETScError::Check(VecSetValues(_vec, 1, &i, &value, INSERT_VALUES));
        }

        void setLocalValue(const PetscInt row, const PetscReal value) {
            PETScError::Check(VecSetValueLocal(_vec, row, value, INSERT_VALUES));
        }

        PetscReal getGlobalValue(const PetscInt i) const {
            PetscReal temp;
            VecGetValues(_vec, 1, &i, &temp);
            return temp;
        }

        PetscReal getGlobalValueWithGhosts(const PetscInt i) const {

            Vec lx;
            VecGhostGetLocalForm(_vec, &lx);

            PetscReal temp;
            VecGetValues(lx, 1, &i, &temp);

            VecGhostRestoreLocalForm(_vec, &lx);
            return temp;
        }

        void addValueGlobal(const PetscInt i, const PetscReal value) {
            VecSetValues(_vec, 1, &i, &value, ADD_VALUES);
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


        void updateGhosts(InsertMode insertmode = INSERT_VALUES) {
            VecGhostUpdateBegin(_vec, insertmode, SCATTER_FORWARD);
            VecGhostUpdateEnd(_vec, insertmode, SCATTER_FORWARD);
        }

        bool hasGhosts() {
            return _nLocalGhosts != 0;
        }


        void describeWithGhosts() {
            Vec lx;
            VecGhostGetLocalForm(_vec, &lx);

            PetscInt n = localSize();
            PetscScalar *array = NULL;
            VecGetArray(lx, &array);

            for (PetscInt i = 0; i < n + _nLocalGhosts; i++) {
                PetscSynchronizedPrintf(_comm, "%D %g\n", i, (double) PetscRealPart(array[i]));
            }

            PetscSynchronizedPrintf(_comm, "\n");
         //   PetscSynchronizedFlush(_comm, PETSC_STDOUT);


            VecRestoreArray(lx, &array);
            VecGhostRestoreLocalForm(_vec, &lx);
        }

        // global size
        PetscInt globalSize() const {
            PetscInt size;
            VecGetSize(_vec, &size);
            return size;
        }

        // local size
        PetscInt localSize() const {
            PetscInt size;
            VecGetLocalSize(_vec, &size);
            return size;
        }

        // global range of the vector
        Range globalRange() const {
            int start, end;
            VecGetOwnershipRange(_vec, &start, &end);
            return Range(start, end);
        }

        // Check if vector is initialized
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

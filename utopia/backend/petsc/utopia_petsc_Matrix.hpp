
#ifndef UTOPIA_UTOPIA_PETSCMATRIX_H
#define UTOPIA_UTOPIA_PETSCMATRIX_H

#include "utopia_petsc_Error.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_Base.hpp"

#include "petscsys.h"
#include "petscmat.h"

#include <memory>

namespace utopia {

    class PETScMatWrapper {
    public:
        PETScMatWrapper(const MPI_Comm comm = PETSC_COMM_WORLD)
                : _comm(comm), owner_(true) {
            MatCreate(_comm, &_mat);
        }

        PETScMatWrapper(Mat &mat, const bool owner = false)
        : _mat(mat), owner_(owner)
        {
            PetscObjectGetComm((PetscObject)_mat,&_comm);
        }

        ~PETScMatWrapper() {
            if(owner_) {
                MatDestroy(&_mat);
            }
        }

        inline Mat &implementation() {
            return _mat;
        }

        inline const Mat &implementation() const {
            return _mat;
        }

        //MAT_DO_NOT_COPY_VALUES or MAT_COPY_VALUES, cause it to copy the numerical values in the matrix MAT_SHARE_NONZERO_PATTERN
        inline void duplicate(PETScMatWrapper &other, MatDuplicateOption opt = MAT_COPY_VALUES) const {
            if(other.owner_) {
                MatDestroy(&other._mat);
            }

            PETScError::Check(MatDuplicate(_mat, opt, &other._mat));
            other._comm = _comm;
        }

        inline void copy(PETScMatWrapper &other, MatStructure opt = DIFFERENT_NONZERO_PATTERN) const {
            //DOES NOT WORK
            assert(false);
            PetscInt rows, cols, grows, gcols;

            MatGetSize(_mat, &grows, &gcols);
            MatGetLocalSize(_mat, &rows, &cols);

            MatSetSizes(other._mat, rows, cols, grows, gcols);

            PETScError::Check(MatCopy(_mat, other._mat, opt));
            other._comm = _comm;
        }

        inline void convert(PETScMatWrapper &other, MatType newtype) {
            //MAT_REUSE_MATRIX is only supported for inplace conversion, otherwise use MAT_INITIAL_MATRIX.
            PETScError::Check(MatConvert(_mat, newtype, MAT_INITIAL_MATRIX, &other._mat));
            other._comm = _comm;
        }

        inline void convert(MatType newtype) {
            //MAT_REUSE_MATRIX is only supported for inplace conversion, otherwise use MAT_INITIAL_MATRIX.
            PETScError::Check(MatConvert(_mat, newtype, MAT_REUSE_MATRIX, &_mat));
        }

        PETScMatWrapper(const PETScMatWrapper &other)
                : _comm(other._comm) {
            //MatStructure str = SAME_NONZERO_PATTERN or DIFFERENT_NONZERO_PATTERN

            PETScError::Check(MatCopy(other._mat, _mat, SAME_NONZERO_PATTERN));
        }

        MPI_Comm &communicator() {
            return _comm;
        }

        void set_owner(const bool owner)
        {
            owner_ = owner;
        }

    private:
        MPI_Comm _comm;
        Mat _mat;
        bool owner_;
    };

    class PETScMatrix {
    public:
        Mat &implementation() {
            return _wrapper->implementation();
        }

        void wrap(Mat &mat)
        {
           _wrapper = std::make_shared<PETScMatWrapper>(mat, false);
        }

        const Mat &implementation() const {
            return _wrapper->implementation();
        }

        PETScMatrix(const MPI_Comm comm = PETSC_COMM_WORLD) {
            using std::make_shared;
            _wrapper = make_shared<PETScMatWrapper>(comm);
        }

        PETScMatrix(const PETScMatrix &other) {
            using std::make_shared;
            _wrapper = make_shared<PETScMatWrapper>();
            other._wrapper->duplicate(*_wrapper);
        }

        PETScMatrix &operator=(const PETScMatrix &other) {
            if(_wrapper == other._wrapper) return *this;

            _wrapper = std::make_shared<PETScMatWrapper>();
            other._wrapper->duplicate(*_wrapper);
            return *this;
        }

        PETScMatrix &operator=(PETScMatrix &&other) {
            if(_wrapper == other._wrapper) return *this;

            this->_wrapper = other._wrapper;
            other._wrapper = nullptr;
            return *this;
        }

        void describe() const {

            MatView(_wrapper->implementation(), PETSC_VIEWER_STDOUT_WORLD);
        }

        MPI_Comm &communicator() {
            return _wrapper->communicator();
        }

        const MPI_Comm &communicator() const {
            return _wrapper->communicator();
        }

    private:
        std::shared_ptr<PETScMatWrapper> _wrapper;
    };
}

#endif //UTOPIA_UTOPIA_PETSCMATRIX_H

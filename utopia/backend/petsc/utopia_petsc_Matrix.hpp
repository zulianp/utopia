
#ifndef UTOPIA_UTOPIA_PETSCMATRIX_H
#define UTOPIA_UTOPIA_PETSCMATRIX_H

#include "utopia_petsc_Error.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_Base.hpp"
#include "utopia_Size.hpp"

#include "petscsys.h"
#include "petscmat.h"

#include <memory>

namespace utopia {
    class PetscMatrixMemory {
    public:
        PetscMatrixMemory(const MPI_Comm comm = PETSC_COMM_WORLD)
        : owner_(true)
        {
            MatCreate(comm, &_mat);
        }

        PetscMatrixMemory(Mat &mat, const bool owner = false)
        : _mat(mat), owner_(owner) { }

        ~PetscMatrixMemory()
        {
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
        inline void duplicate(PetscMatrixMemory &other, MatDuplicateOption opt = MAT_COPY_VALUES) const {
            if(other.owner_) {
                MatDestroy(&other._mat);
            }

            PetscError::Check(MatDuplicate(_mat, opt, &other._mat));
        }

        inline void convert(PetscMatrixMemory &other, MatType newtype) {
            //MAT_REUSE_MATRIX is only supported for inplace conversion, otherwise use MAT_INITIAL_MATRIX.
            PetscError::Check(MatConvert(_mat, newtype, MAT_INITIAL_MATRIX, &other._mat));
        }

        inline void convert(MatType newtype) {
            //MAT_REUSE_MATRIX is only supported for inplace conversion, otherwise use MAT_INITIAL_MATRIX.
            PetscError::Check(MatConvert(_mat, newtype, MAT_REUSE_MATRIX, &_mat));
        }

        PetscMatrixMemory(const PetscMatrixMemory &other)
        : owner_(true)
        {
            PetscError::Check(MatCopy(other._mat, _mat, SAME_NONZERO_PATTERN));
        }

       inline MPI_Comm communicator() const {
           MPI_Comm comm = PetscObjectComm((PetscObject) implementation());
           assert(comm != MPI_COMM_NULL);
           return comm;
       }

        void set_owner(const bool owner)
        {
            owner_ = owner;
        }

     

    private:
        Mat _mat;
        bool owner_;
    };

    class PetscMatrix {
    public:
        Mat &implementation() {
            return _wrapper->implementation();
        }

        void wrap(Mat &mat)
        {
           _wrapper = std::make_shared<PetscMatrixMemory>(mat, false);
        }

        const Mat &implementation() const {
            return _wrapper->implementation();
        }

        PetscMatrix(const MPI_Comm comm = PETSC_COMM_WORLD) {
            using std::make_shared;
            _wrapper = make_shared<PetscMatrixMemory>(comm);
        }

        PetscMatrix(const PetscMatrix &other) {
            using std::make_shared;
            _wrapper = make_shared<PetscMatrixMemory>();
            other._wrapper->duplicate(*_wrapper);
        }

        PetscMatrix &operator=(const PetscMatrix &other) {
            if(_wrapper == other._wrapper) return *this;

            _wrapper = std::make_shared<PetscMatrixMemory>();
            other._wrapper->duplicate(*_wrapper);
            return *this;
        }

        PetscMatrix &operator=(PetscMatrix &&other) {
            if(_wrapper == other._wrapper) return *this;

            this->_wrapper = other._wrapper;
            other._wrapper = nullptr;
            return *this;
        }

        inline void describe() const {
            MatView(_wrapper->implementation(), PETSC_VIEWER_STDOUT_(communicator()));
        }

        inline MPI_Comm communicator() const {
            return _wrapper->communicator();
        }

        inline Size size() const
        {
            PetscInt rows, cols;
            MatGetSize(implementation(), &rows, &cols);
            return {rows, cols};
        }

        inline Size local_size() const
        {
            PetscInt rows, cols;
            MatGetLocalSize(implementation(), &rows, &cols);
            return {rows, cols};
        }

    private:
        std::shared_ptr<PetscMatrixMemory> _wrapper;
    };
}

#endif //UTOPIA_UTOPIA_PETSCMATRIX_H

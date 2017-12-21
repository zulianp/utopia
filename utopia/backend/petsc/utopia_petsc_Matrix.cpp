#include "utopia_petsc_Matrix.hpp"

// /*
//     matrix impl in the hpp file
//     real    0m11.091s
//     user    0m6.140s
//     sys 0m4.893s


//     this implementation
//     real    0m35.874s
//     user    0m30.461s
//     sys 0m5.210s
// */


// #ifndef UTOPIA_UTOPIA_PETSCMATRIX_H
// #define UTOPIA_UTOPIA_PETSCMATRIX_H

// #include "utopia_petsc_Error.hpp"
// #include "utopia_petsc_Vector.hpp"
// #include "utopia_Base.hpp"

// #include "petscsys.h"
// #include "petscmat.h"

// #include <memory>

// namespace utopia {

//     class PETScMatrix {
//     public:
//         inline PETScMatrix(const MPI_Comm comm = PETSC_COMM_WORLD)
//         : owner_(true)
//         {
//             MatCreate(comm, &mat_);
//         }

//         PETScMatrix(const PETScMatrix &other)
//         : owner_(true)
//         {
//             PETScError::Check(MatDuplicate(other.mat_, MAT_COPY_VALUES, &mat_));
//         }

//         inline PETScMatrix(PETScMatrix &&other) 
//         : mat_(other.mat_), owner_(other.owner_)
//         {
//             other.owner_ = false;
//             other.mat_ = nullptr;
//         }

//         inline PETScMatrix(Mat &mat, const bool owner = false)
//         : mat_(mat), owner_(owner) { }

//         inline PETScMatrix &operator=(const PETScMatrix &other) {
//             if(this->mat_ == other.mat_) return *this;
//             other.duplicate(*this);
//             return *this;
//         }

//         inline PETScMatrix &operator=(PETScMatrix &&other) {
//             if(this->mat_ == other.mat_) return *this;

//             mat_ = other.mat_;
//             owner_ = other.owner_;

//             other.owner_ = false;
//             other.mat_ = nullptr;
//             return *this;
//         }

//         inline PETScMatrix &operator=(const PETScMatrix &&other) {
//             if(this->mat_ == other.mat_) return *this;

//             assert(false);
//             return *this;
//         }

//         ~PETScMatrix()
//         {
//             destroy();
//         }

//         inline void destroy()
//         {
//             if(owner_) {
//                 MatDestroy(&mat_);
//             }
//         }

//         void wrap(Mat &mat)
//         {
//             destroy();

//             mat_ = mat; 
//             owner_ = false;
//         }

//         inline Mat &implementation() {
//             return mat_;
//         }

//         inline const Mat &implementation() const {
//             return mat_;
//         }

//         inline void convert(PETScMatrix &other, MatType newtype) {
//             PETScError::Check(MatConvert(mat_, newtype, MAT_INITIAL_MATRIX, &other.mat_));
//         }

//         inline void convert(MatType newtype) {
//             PETScError::Check(MatConvert(mat_, newtype, MAT_REUSE_MATRIX, &mat_));
//         }

//         inline MPI_Comm communicator() const {
//             MPI_Comm comm = PetscObjectComm((PetscObject) implementation());
//             assert(comm != MPI_COMM_NULL);
//             return comm;
//         }

//         void set_owner(const bool owner)
//         {
//             owner_ = owner;
//         }

//         inline void describe() const {
//             MatView(implementation(), PETSC_VIEWER_STDOUT_WORLD);
//         }

//     private:
//         inline void duplicate(PETScMatrix &other, MatDuplicateOption opt = MAT_COPY_VALUES) const {
//             other.destroy();
//             PETScError::Check(MatDuplicate(mat_, opt, &other.mat_));
//         }

//         Mat mat_;
//         bool owner_;
//     };
// }

// #endif //UTOPIA_UTOPIA_PETSCMATRIX_H



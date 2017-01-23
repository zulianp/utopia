
#ifndef UTOPIA_UTOPIA_PETSCMATRIX_H
#define UTOPIA_UTOPIA_PETSCMATRIX_H


#include "utopia_PETScError.hpp"
#include "utopia_PETScVector.hpp"
#include "utopia_Base.hpp"

#include "petscsys.h"
#include "petscmat.h"

#include <memory>

namespace utopia {


    class PETScMatrixMeta {
    public:
        void setLocalSize(const PetscInt localRows, const PetscInt localColumns) {
            _localRows = localRows;
            _localColumns = localColumns;
        }

        void setGlobalSize(const PetscInt globalRows, const PetscInt globalColumns) {
            _globalRows = globalRows;
            _globalColumns = globalColumns;
        }

        void setInsertMode(InsertMode insertMode) {
            _insertMode = insertMode;
        }

        std::vector<PetscInt> &nnzXRow() {
            return _nnzXrow;
        }

        const std::vector<PetscInt> &nnzXRow() const {
            return _nnzXrow;
        }

        inline PetscInt getGlobalRows() const {
            return _globalRows;
        }

        inline PetscInt getGlobalColumns() const {
            return _globalColumns;
        }


        inline PetscInt getLocalRows() const {
            return _localRows;
        }

        inline PetscInt getLocalColumns() const {
            return _localColumns;
        }

        inline InsertMode getInsertMode() const {
            return _insertMode;
        }

        PETScMatrixMeta()
                : _globalRows(PETSC_DETERMINE), _globalColumns(PETSC_DETERMINE),
                  _localRows(PETSC_DECIDE), _localColumns(PETSC_DECIDE),
                  _insertMode(INSERT_VALUES) { }

        void describe(std::ostream &os, MPI_Comm comm = PETSC_COMM_WORLD) {

            int rank;
            int size;

            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &size);

            MPI_Barrier(comm);

            for (int r = 0; r < size; ++r) {
                if (r == rank) {
                    os << "On rank " << rank << "\n";
                    os << "Global size " << _globalRows << " x " << _globalColumns << "\n";
                    os << "Local size " << _localRows << " x " << _localColumns << "\n";
                    os << "Nnz x row\n";
                    disp(_nnzXrow.begin(), _nnzXrow.end(), os);

                }

                MPI_Barrier(comm);
            }
        }

    private:
        PetscInt _globalRows, _globalColumns;
        PetscInt _localRows, _localColumns;
        std::vector<PetscInt> _nnzXrow;
        InsertMode _insertMode;

    };

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

        bool scaleByDiagonalMatrixFromLeft(PETScVector &vector)
        {
            return PETScError::Check( MatDiagonalScale(_wrapper->implementation(), vector.implementation(), NULL) );
        }

        bool mul(PETScMatrix &rhs, PETScMatrix &result) {
            return PETScError::Check(
                    MatMatMult(_wrapper->implementation(), rhs._wrapper->implementation(), MAT_INITIAL_MATRIX,
                               PETSC_DEFAULT, &result._wrapper->implementation()));
        }

        bool mul(PETScVector &rhs, PETScVector &result) {
            if (!result.hasGhosts()) {

                PetscInt globalRows, globalColumns;
                MatGetSize(_wrapper->implementation(), &globalRows, &globalColumns);
                PetscInt localRows, localColumns;
                MatGetLocalSize(_wrapper->implementation(), &localRows, &localColumns);

                result.initialize(localRows, globalRows);
            }

            MatMult(_wrapper->implementation(), rhs.implementation(), result.implementation());
            return true;
        }

        inline Range ownershipRange() const {
            PetscInt globalBegin, globalEnd;
            MatGetOwnershipRange(_wrapper->implementation(), &globalBegin, &globalEnd);
            return Range(globalBegin, globalEnd);
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

        bool get(const std::vector<PetscInt> &rowIndex, const std::vector<PetscInt> &colIndex,
                 std::vector<PetscReal> &values) {
            values.resize(rowIndex.size());
            return PETScError::Check(
                    MatGetValues(_wrapper->implementation(), rowIndex.size(), &rowIndex[0], colIndex.size(),
                                 &colIndex[0], &values[0]));
        }

        inline bool set(const std::vector<PetscInt> &rowIndex,
                        const std::vector<PetscInt> &colIndex,
                        const std::vector<PetscReal> &values) {
            return insert(rowIndex, colIndex, values, INSERT_VALUES);
        }

        inline bool add(const std::vector<PetscInt> &rowIndex,
                        const std::vector<PetscInt> &colIndex,
                        const std::vector<PetscReal> &values) {
            return insert(rowIndex, colIndex, values, ADD_VALUES);
        }

        void describe() const {

            MatView(_wrapper->implementation(), PETSC_VIEWER_STDOUT_WORLD);
        }

        bool diag(PETScVector &vec) {
            using std::max;
            PetscInt globalRows, globalColumns;
            MatGetSize(_wrapper->implementation(), &globalRows, &globalColumns);
            PetscInt localRows, localColumns;
            MatGetLocalSize(_wrapper->implementation(), &localRows, &localColumns);

            PetscInt lenGlobal = max(globalRows, globalColumns);
            PetscInt lenLocal = max(localRows, localColumns);
            vec.resize(lenLocal, lenGlobal);
            return PETScError::Check(MatGetDiagonal(_wrapper->implementation(), vec.implementation()));
        }

        // bool setDiag(PETScVector &vec)
        // {
        // 	return PETScError::Check( MatSetDiagonal(_wrapper->implementation(), vec.implementation() ) );
        // }


        void save(const std::string &path) {
            PetscViewer viewer;

            PetscViewerCreate(_wrapper->communicator(), &viewer);
            PetscViewerSetType(viewer, PETSCVIEWERASCII);
            PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
            PetscViewerFileSetName(viewer, path.c_str());
            PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
            MatView(_wrapper->implementation(), viewer);
            PetscViewerDestroy(&viewer);
        }

        void setName(const std::string &name) {
            PetscObjectSetName((PetscObject) _wrapper->implementation(), name.c_str());
        }

//        inline Structure localStructure() const
//        {
//            PetscInt localRows, localColumns;
//            MatGetLocalSize(_wrapper->implementation(), &localRows, &localColumns);
//            return Structure(localRows, localColumns);
//        }
//
//        inline Structure globalStructure() const
//        {
//            PetscInt globalRows, globalColumns;
//            MatGetSize(_wrapper->implementation(), &globalRows, &globalColumns);
//            return Structure(globalRows, globalColumns);
//        }

        MPI_Comm &communicator() {
            return _wrapper->communicator();
        }

        const MPI_Comm &communicator() const {
            return _wrapper->communicator();
        }

//    private:
        bool finalize() {
            bool ok = PETScError::Check(MatAssemblyBegin(_wrapper->implementation(), MAT_FINAL_ASSEMBLY));
            ok &= PETScError::Check(MatAssemblyEnd(_wrapper->implementation(), MAT_FINAL_ASSEMBLY));
            return ok;
        }

        inline bool insert(const std::vector<PetscInt> &rowIndex,
                           const std::vector<PetscInt> &colIndex,
                           const std::vector<PetscReal> &values,
                           InsertMode addv) {
            // std::cout << "Inserting:\n r:" << rowIndex << "\n c:" << colIndex << "\n v:" << values << std::endl;
            assert(rowIndex.size() * colIndex.size() == values.size());
            assert(rowIndex.size() >= 1);
            assert(colIndex.size() >= 1);

            return PETScError::Check(MatSetValues(_wrapper->implementation(), rowIndex.size(), &rowIndex[0],
                                                  colIndex.size(), &colIndex[0], &values[0],
                                                  addv));
        }

        inline bool initialize(const PetscInt localRows, const PetscInt localColumns,
                               const PetscInt /*nDiags */,
                               std::vector<PetscInt> &nnzXrow,
                               std::vector<PetscInt> &rowIndex,
                               std::vector<PetscInt> &colIndex,
                               std::vector<PetscReal> &values,
                               InsertMode addv = INSERT_VALUES) {
            // MatMPIAIJSetPreallocation(_wrapper->implementation(), values.size(), NULL, 5,NULL);
            initLocal(localRows, localColumns);
            MatSetType(_wrapper->implementation(), MATMPIAIJ);
            // MatSetUp(_wrapper->implementation());//, PETSC_DEFAULT, NULL, values.size(), NULL);/

            // MatCreateAIJ(_wrapper->communicator(), localRows, localColumns, PETSC_DETERMINE, PETSC_DETERMINE, nDiags, NULL, values.size(), NULL, &_wrapper->implementation());
            MatSeqAIJSetPreallocation(_wrapper->implementation(), PETSC_DEFAULT, &nnzXrow[0]);
            
            // MatMPIAIJSetPreallocation(_wrapper->implementation(), nDiags, NULL, values.size()-nDiags, NULL);
            MatMPIAIJSetPreallocation(_wrapper->implementation(), PETSC_DEFAULT, NULL, PETSC_DEFAULT, NULL);

            // initLocal(localRows, localColumns);
            bool ok = PETScError::Check(MatAssemblyBegin(_wrapper->implementation(), MAT_FLUSH_ASSEMBLY));
            ok &= PETScError::Check(
                    MatSetValues(_wrapper->implementation(), rowIndex.size(), &rowIndex[0], colIndex.size(),
                                 &colIndex[0], &values[0], addv));
            ok &= PETScError::Check(MatAssemblyEnd(_wrapper->implementation(), MAT_FLUSH_ASSEMBLY));

            return finalize();
        }


        void initGlobal(const PetscInt globalRows, const PetscInt globalColumns) {
            MatSetSizes(_wrapper->implementation(), PETSC_DECIDE, PETSC_DECIDE, globalRows, globalColumns);
        }

        void initLocal(const PetscInt localRows, const PetscInt localColumns) {
            MatSetSizes(_wrapper->implementation(), localRows, localColumns, PETSC_DETERMINE, PETSC_DETERMINE);
        }

    private:
        std::shared_ptr<PETScMatWrapper> _wrapper;

    };
}

#endif //UTOPIA_UTOPIA_PETSCMATRIX_H

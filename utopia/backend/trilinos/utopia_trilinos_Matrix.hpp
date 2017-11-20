
#ifndef UTOPIA_TPETRAMATRIX_H
#define UTOPIA_TPETRAMATRIX_H

#include "utopia_trilinos.hpp"
#include "utopia_trilinos_Matrix.hpp"

#include "utopia_Base.hpp"

#include <memory>

namespace utopia {
    typedef Tpetra::CrsMatrix<>                       crs_matrix_type;
    typedef Teuchos::RCP<crs_matrix_type>             rcp_crs_matrix_type;
    typedef Teuchos::RCP<const Teuchos::Comm<int> >   rcp_comm_type;
    typedef Tpetra::Map<>                             map_type;
    typedef Teuchos::RCP< const map_type >            rcp_map_type;
    typedef Tpetra::Vector<>::local_ordinal_type      local_ordinal_type;
    typedef Tpetra::Vector<>::global_ordinal_type     global_ordinal_type;
    typedef Tpetra::Vector<>::scalar_type             scalar_type;
    class TpetraMatrixMeta {
    public:

        void setLocalSize(const int localRows, const int localColumns) {
            _localRows = localRows;
            _localColumns = localColumns;
        }

        void setGlobalSize(const int globalRows, const int globalColumns) {
            _globalRows = globalRows;
            _globalColumns = globalColumns;
        }

        std::vector<int> &nnzXRow() {
            return _nnzXrow;
        }

        const std::vector<int> &nnzXRow() const {
            return _nnzXrow;
        }

        inline int getGlobalRows() const {
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

/*        TpetraMatrixMeta()
                : _globalRows(TPETRA_DETERMINE), _globalColumns(TPETRA_DETERMINE),
                  _localRows(TPETRA_DECIDE), _localColumns(TPETRA_DECIDE),
                  _insertMode(INSERT_VALUES) { }
*/



        void describe(std::ostream &os, MPI_Comm comm ) {

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
        };

    private:

        int _globalRows; 
        int _globalColumns;
        int _localRows; 
        int _localColumns;
        std::vector<int> _nnzXrow;

    };//TpetraMatrixMeta

    class TpetraMatWrapper {
    public:
        TpetraMatWrapper(const rcp_comm_type comm )
                : _comm(comm), _owner(true) {

        Teuchos::RCP<const map_type> _map = Teuchos::rcp (new map_type (10, 0, _comm));

        Teuchos::RCP<crs_matrix_type> _mat (new crs_matrix_type (_map, 0)); 
        }

        TpetraMatWrapper()
            {
            _comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
            _map = Teuchos::rcp (new Tpetra::Map<> (_comm->getSize (), 0, _comm));
            _mat.reset(new crs_matrix_type (_map,0));
            _owner = true;
           }

        //copy constructor based on Trilinos matrix
        TpetraMatWrapper(crs_matrix_type  matrix)
            {       
            _map = matrix.getMap();
            _comm = _map->getComm ();
         //   _mat.reset (new crs_matrix_type (matrix,Teuchos::View));  //TODO copy constructor
            _owner = true;
            }

        TpetraMatWrapper(rcp_map_type  map)
            {
            _comm = map->getComm ();
            _map = map;
            _mat.reset ( new crs_matrix_type (_map,0));
            _owner = true;
        }

      ~TpetraMatWrapper() {
            if(_owner) {
                //MatDestroy(&_mat);
            }
        }

        inline rcp_crs_matrix_type &implementation() {
            return _mat;
        }

        inline const rcp_crs_matrix_type &implementation() const {
            return _mat;
        }

        inline bool insertGlobalValues(global_ordinal_type gblRow, Teuchos::tuple<global_ordinal_type> (gblRow, gblRow + 1), Teuchos::tuple<scalar_type> (two, negOne)){
        _mat->insertGlobalValues(gblRow, Teuchos::tuple<global_ordinal_type> (gblRow, gblRow + 1), Teuchos::tuple<scalar_type> (two, negOne) );
        return true;
        }

        inline bool fillComplete(){
        _mat->fillComplete();
        return true;
        }

        //MAT_DO_NOT_COPY_VALUES or MAT_COPY_VALUES, cause it to copy the numerical values in the matrix MAT_SHARE_NONZERO_PATTERN
  /*      inline void duplicate(PETScMatWrapper &other, MatDuplicateOption opt = MAT_COPY_VALUES) const {
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
*/
        rcp_comm_type &communicator() {
            return _comm;
        }

        void set_owner(const bool owner)
        {
            _owner = owner;
        }

    private:
        rcp_comm_type _comm;
        rcp_crs_matrix_type  _mat;
        Teuchos::RCP<const map_type> _map;
        bool _owner;
    };  //TpetraMatWrapper

    class TpetraMatrix {
    public:
        TpetraMatrix()
            {
           _rcp_wrapper.reset(new TpetraMatWrapper() );
           }

        //copy constructor based on Trilinos matrix
      /*  TpetraMatrix(matrix_type  matrix)
            {  
            _rcp_wrapper.reset(new TpetraMatWrapper(matrix) );    //TODO
            }*/

        TpetraMatrix(rcp_map_type  map)
            {
            _rcp_wrapper.reset(new TpetraMatWrapper(map) );
        }

        rcp_crs_matrix_type  &implementation() {
            return _rcp_wrapper->implementation();
        }

        void wrap(rcp_crs_matrix_type &mat)
        {
           //_rcp_wrapper = std::make_shared<TpetraMatWrapper>(mat, false); //TODO overloading of =
        }

        const rcp_crs_matrix_type &implementation() const {
            return _rcp_wrapper->implementation();
        }

 /*       bool scaleByDiagonalMatrixFromLeft(PETScVector &vector)
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

        TpetraMatrix(const rcp_comm_type comm ) {
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
            // _wrapper = std::make_shared<PETScMatWrapper>();
            // other._wrapper->duplicate(*_wrapper);
            return *this;
        }

        bool get(const std::vector<PetscInt> &rowIndex, const std::vector<PetscInt> &colIndex,
                 std::vector<PetscReal> &values) {
            values.resize(rowIndex.size());
            return PETScError::Check(
                    MatGetValues(_wrapper->implementation(), rowIndex.size(), &rowIndex[0], colIndex.size(),
                                 &colIndex[0], &values[0]));
        }
*/
        inline bool set(const std::vector<int> &rowIndex,
                        const std::vector<int> &colIndex,
                        const std::vector<double> &values) {
            return insert(rowIndex, colIndex, values);//, INSERT_VALUES);
        }

        inline bool add(const std::vector<int> &rowIndex,
                        const std::vector<int> &colIndex,
                        const std::vector<double> &values) {
            return insert(rowIndex, colIndex, values);//, ADD_VALUES);
        }
/*
        void describe() const {

            MatView(_wrapper->implementation(), PETSC_VIEWER_STDOUT_WORLD);
        }

        bool diag(PETScVector &vec) {
            using std::max;
            int globalRows, globalColumns;
            MatGetSize(_wrapper->implementation(), &globalRows, &globalColumns);
            int localRows, localColumns;
            MatGetLocalSize(_wrapper->implementation(), &localRows, &localColumns);

            int lenGlobal = max(globalRows, globalColumns);
            int lenLocal = max(localRows, localColumns);
            vec.resize(lenLocal, lenGlobal);
            return PETScError::Check(MatGetDiagonal(_wrapper->implementation(), vec.implementation()));
        }
*/
        // bool setDiag(PETScVector &vec)
        // {
        // 	return PETScError::Check( MatSetDiagonal(_wrapper->implementation(), vec.implementation() ) );
        // }


        rcp_comm_type &communicator() {
            return _rcp_wrapper->communicator();
        }

        const rcp_comm_type &communicator() const {
            return _rcp_wrapper->communicator();
        }

//    private:
        bool finalize() {//fiilComplete
            /*bool ok = PETScError::Check(MatAssemblyBegin(_wrapper->implementation(), MAT_FINAL_ASSEMBLY));
            ok &= PETScError::Check(MatAssemblyEnd(_wrapper->implementation(), MAT_FINAL_ASSEMBLY));
            return ok;*/
        if _rcp_wrapper->fillComplete()
        return true;
        }

        inline bool insert(const std::vector<int> &rowIndex,
                           const std::vector<int> &colIndex,
                           const std::vector<double> &values){
                           //InsertMode addv) {
            std::cout << "Inserting:\n r:" << rowIndex << "\n c:" << colIndex << "\n v:" << values << std::endl;
            assert(rowIndex.size() * colIndex.size() == values.size());
            assert(rowIndex.size() >= 1);
            assert(colIndex.size() >= 1);
return false;
//     return PETScError::Check(MatSetValues(_wrapper->implementation(), rowIndex.size(), &rowIndex[0],colIndex.size(), &colIndex[0], &values[0], addv));
        }



    private:
        Teuchos::RCP<TpetraMatWrapper> _rcp_wrapper;

    }; //TpetraMatrix
}

#endif //UTOPIA_TPETRAMATRIX_H

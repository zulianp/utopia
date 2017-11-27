
#ifndef UTOPIA_TPETRAVECTOR_H
#define UTOPIA_TPETRAVECTOR_H

#include "utopia_trilinos.hpp"
//#include "utopia_trilinos_Types.hpp"

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"

#include <memory>

namespace utopia {

    typedef Tpetra::Map<>                             map_type;
    typedef Tpetra::Vector<>                          vector_type;

    typedef vector_type::global_ordinal_type            global_ordinal_type;

    typedef Teuchos::RCP<const Teuchos::Comm<int> >     rcp_comm_type;
    typedef Teuchos::RCP<const map_type>                rcp_map_type;
    typedef Teuchos::RCP<const vector_type>             rcp_vector_type;
    typedef Tpetra::CrsMatrix<>                         crs_matrix_type;
    typedef Teuchos::RCP<crs_matrix_type>               rcp_crs_matrix_type;

        typedef Tpetra::Vector<>::scalar_type         scalar_type;
        typedef Tpetra::Vector<>::local_ordinal_type  local_ordinal_type;
        typedef Tpetra::Vector<>::mag_type            magnitude_type;


	
    class TpetraVector {

    public:
        TpetraVector()
            {
            // FIXME global size to size of the comm and index base to zero
            _comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
            _contigMap = Teuchos::rcp (new Tpetra::Map<> (_comm->getSize (), 0, _comm));
            _vec.reset(new vector_type (_contigMap));
            _initialized = false; 
            _nLocalGhosts = 0;
            }

        //copy constructor based on Trilinos vector
        TpetraVector(vector_type  vector)
            {          
            _contigMap = vector.getMap();
            _comm = _contigMap->getComm ();
            _vec.reset (new vector_type (vector, Teuchos::View));
            _initialized = false;
            _nLocalGhosts = 0;
            }

        TpetraVector(rcp_map_type  map)
            {
            // FIXME global size to size of the comm and index base to zero
            _comm = map->getComm ();
            _contigMap = map;
            _vec.reset ( new vector_type (_contigMap));
            _initialized = false; 
            _nLocalGhosts = 0;
            }

/*        TpetraVector(const utopia::TpetraVector&  vector )
            {
            _contigMap = vector.getMap();
            _comm = vector.getComm ();
            _vec.reset (new vector_type (vector, Teuchos::View));
            _initialized = vector.isInitialized();
            _nLocalGhosts = vector.getLocalGhosts();
            }*/


        ~TpetraVector() {
            destroy();
        }

/*        TpetraVector(rcp_map_type  map, bool init=false, int ghost=0)
		    : _comm(map->getComm ()), _contigMap(map), _initialized(init), _nLocalGhosts(ghost), _vec  (new const vector_type (map)){
            // FIXME global size to size of the comm and index base to zero
           // _vec(_contigMap);
          //_vec = new const vector_type (map);
//          vector_type _vec (map);
        }*/



        TpetraVector(const TpetraVector &other) {
            // Copy communicator
            _comm = other._comm;
            // Destory current vector

            // Copy other
           // PETScError::Check(VecDuplicate(other._vec, &_vec));
           // PETScError::Check(VecCopy(other._vec, _vec));
            // Copy fields
            // FIXME: Need to copy this fields??
            _initialized = other._initialized;
            _nLocalGhosts = other._nLocalGhosts;
            _ghosts = other._ghosts;
        }

        // get communicator
        rcp_comm_type &communicator() {
            return _comm;
        }

        const rcp_comm_type &communicator() const {
            return _comm;
        }

        // set communicator for MPI
        void setCommunicator(const rcp_comm_type comm) {
            _comm = comm;
        }

        // assign operator
        TpetraVector &operator=(const TpetraVector &other) {
            if(this == &other) return *this;
            // TODO
            //_wrapper = std::make_shared<PETScVectorWrapper>();
            //other._wrapper->copy(*_wrapper);
            destroy();
         //   PETScError::Check(VecDuplicate(other._vec, &_vec));
         //   PETScError::Check(VecCopy(other._vec, _vec));
            _comm = other._comm;

            return *this;
        }

        TpetraVector &operator=(TpetraVector &&other) {
            if(this == &other) return *this;
            // TODO
            //_wrapper = std::make_shared<PETScVectorWrapper>();
            //other._wrapper->copy(*_wrapper);
            destroy();
            // PETScError::Check(VecDuplicate(other._vec, &_vec));
            // PETScError::Check(VecCopy(other._vec, _vec));
            _comm = other._comm;
            //_vec= other._vec;
            //other._vec = nullptr; //TODO

            return *this;
        }

        // initialize vector
        bool initialize(const int nLocal, const int nGlobal) {
            // FIXME Error check?
         //   PETScError::Check(VecCreateMPI(_comm, nLocal, nGlobal, &_vec));
            return true;
        }

        const std::vector<int> &ghosts() const {
            return _ghosts;
        }

        bool initializeWithGhosts(const int nLocal, const int nGlobal, const std::vector<int> &ghosts) {
            // FIXME Error check?
        //    PETScError::Check(VecCreate(PETSC_COMM_WORLD, &_vec));
        //    PETScError::Check(VecSetType(_vec, VECMPI));
         //   PETScError::Check(VecSetSizes(_vec, nLocal, nGlobal));

            if (!ghosts.empty()) {
         //       VecMPISetGhost(_vec, ghosts.size(), &ghosts[0]);
                _nLocalGhosts = ghosts.size();
                _ghosts = ghosts;
                //updateGhosts();
            }

            return true;
        }

        // destroy vector
        void destroy() {
            //VecDestroy(&_vec);
        }

        void resize(const int nLocal, const int nGlobal ) {
            destroy();
            // FIXME Error check?
          //  PETScError::Check(VecCreateMPI(_comm, nLocal, nGlobal, &_vec));
            _initialized = true;
        }

        // assembly vector
        void assemblyBegin() {
        //    PETScError::Check(VecAssemblyBegin(_vec));
        }

        void assemblyEnd() {
        //    PETScError::Check(VecAssemblyEnd(_vec));
        }

        void finalize() {
            assemblyBegin();
            assemblyEnd();
        }

        void setGlobalValue(const int i, const double value) {
        //    PETScError::Check(VecSetValues(_vec, 1, &i, &value, INSERT_VALUES));
        }

        void setLocalValue(const int row, const double value) {
        //    PETScError::Check(VecSetValueLocal(_vec, row, value, INSERT_VALUES));
        }

        double getGlobalValue(const int i) const {
            double temp;
           // VecGetValues(_vec, 1, &i, &temp);
            return temp;
        }

        double getGlobalValueWithGhosts(const int i) const {

            //vector_type lx;
          //  VecGhostGetLocalForm(_vec, &lx);

            double temp;
          //  VecGetValues(lx, 1, &i, &temp);

          //  VecGhostRestoreLocalForm(_vec, &lx);
            return temp;
        }

     /*   void addValueGlobal(const int i, const double value) {
            VecSetValues(_vec, 1, &i, &value, ADD_VALUES);
        }*/

     /*   vector_type &implementation() {
            return _vec;
        }

        const vector_type &implementation() const {
            return _vec;
        }*/

       /* void describe() const {
          //TODO kokkos view
          //  VecView(_vec, PETSC_VIEWER_STDOUT_WORLD);
        }*/


     /*   void updateGhosts(InsertMode insertmode = INSERT_VALUES) {
            VecGhostUpdateBegin(_vec, insertmode, SCATTER_FORWARD);
            VecGhostUpdateEnd(_vec, insertmode, SCATTER_FORWARD);
        }*/

        bool hasGhosts() {
            return _nLocalGhosts != 0;
        }

        int getLocalGhosts() {
            return _nLocalGhosts;
        }

        rcp_map_type getMap() {
            return _contigMap;
        }

        rcp_comm_type getComm() {
            return _comm;
        }

    /*    void describeWithGhosts() {
            vector_type lx;
            VecGhostGetLocalForm(_vec, &lx);

            int n = localSize();
            double *array = NULL;
            VecGetArray(lx, &array);

            for (int i = 0; i < n + _nLocalGhosts; i++) {
            //    PetscSynchronizedPrintf(_comm, "%D %g\n", i, (double) doublePart(array[i]));
            }

          //  PetscSynchronizedPrintf(_comm, "\n");


            VecRestoreArray(lx, &array);
            VecGhostRestoreLocalForm(_vec, &lx);
        }*/

        // global size
    /*    int globalSize() const {
            int  size;
       //     VecGetSize(_vec, &size);
            return size;
        } */

        // local size
     /*   int localSize() const {
            int size;
        //get local size
           // VecGetLocalSize(_vec, &size);
            return size;
        }*/

        // global range of the vector
        Range globalRange() const {
            int start, end;
            //VecGetOwnershipRange(_vec, &start, &end);
            return Range(start, end);
        }

        // Check if vector is initialized
        inline bool isInitialized() const {
            return _initialized;
        }


    // private data
    private:
        rcp_comm_type  _comm;
        rcp_map_type   _contigMap;
        rcp_vector_type    _vec;

        bool _initialized;
        int _nLocalGhosts;
        std::vector<int> _ghosts;
    };
}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H

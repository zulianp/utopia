#ifndef UTOPIA_TPETRAMATRIX_H
#define UTOPIA_TPETRAMATRIX_H

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"
#include "utopia_Size.hpp"
#include "utopia_Logger.hpp"

#include "utopia_Tpetra_Vector.hpp"

#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_Map_decl.hpp>


#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>

#include <iostream>
#include <memory>

namespace utopia {
    //template<class NodeType>
    class TpetraMatrix {
    public:

        /////////////////////////////////////////////////////////////
        // typedef definitions
        /////////////////////////////////////////////////////////////

        //types of Operators
        typedef Tpetra::Operator<>::scalar_type SC;
        typedef Tpetra::Operator<SC>::local_ordinal_type LO;
        typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;

        //types of Kokkos Parallel Nodes
        typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;
#ifdef KOKKOS_CUDA
        typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
        typedef cuda_node NT;
#elif defined KOKKOS_OPENMP
        typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
        typedef Kokkos::Compat::KokkosThreadsWrapperNode thread_node;
        typedef openmp_node NT;
#else
        typedef serial_node NT;
#endif
        //types of Trilinos Objects
        typedef Tpetra::CrsMatrix<SC, LO, GO, NT>             crs_mat_type;
        typedef Teuchos::RCP<crs_mat_type>                    rcp_crs_mat_type;
        typedef Teuchos::RCP<const Teuchos::Comm<int> >       rcp_comm_type;
        typedef Tpetra::Map<LO, GO, NT>                       map_type;
        typedef Teuchos::RCP<const map_type>                  rcp_map_type;
        typedef Tpetra::Vector<SC, LO, GO, NT>::scalar_type   Scalar;

        /////////////////////////////////////////////////////////////
        //Constructors
        /////////////////////////////////////////////////////////////

        //Default Constructor
        TpetraMatrix() : owner_(true) {}
        TpetraMatrix(int Ndofs, int maxNumEntries) : owner_(true) {
        int indexBase = 0;
        rcp_comm_type Comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
        rcp_map_type Map = Teuchos::rcp(new map_type(Ndofs, indexBase, Comm));
        mat_.reset(new crs_mat_type(Map, maxNumEntries, Tpetra::StaticProfile));
        }

        //Explicit Constructors
        TpetraMatrix(rcp_map_type Map, int maxNumEntries) : owner_(true) {
        mat_.reset(new crs_mat_type(Map, maxNumEntries, Tpetra::StaticProfile));
        }

        //deep copy
        //     template <class Node2>
        // rcp_crs_mat_type  clone (
        //   const Teuchos::RCP<Node2> & node2,
        //   const Teuchos::RCP<Teuchos::ParameterList> & params = Teuchos::null)
        // {
        //   return mat_->clone(node2, params);
        // }

        TpetraMatrix(const TpetraMatrix &other)
        : owner_(true)
        {
            if(!other.is_null()) {
                mat_ = other.mat_->clone(other.mat_->getNode());
            }
        }

        TpetraMatrix(TpetraMatrix &&other)
        : mat_(std::move(other.mat_)), owner_(std::move(other.owner_))
        {}

        /////////////////////////////////////////////////////////////
        //Destructor
        /////////////////////////////////////////////////////////////

         ~TpetraMatrix()
         {}

        /////////////////////////////////////////////////////////////
        //Overloading Operators
        /////////////////////////////////////////////////////////////

        TpetraMatrix &operator=(const TpetraMatrix &other)
        {
            if(this == &other) return *this;

            if(other.is_null()) {
                mat_.reset();
                owner_ = true;
                return *this;
            }

            mat_ = other.mat_->clone(other.mat_->getNode());
            owner_ = true;
            return *this;
        }

        TpetraMatrix &operator=(TpetraMatrix &&other)
        {
            if(this == &other) return *this;

            if(other.is_null()) {
                mat_.reset();
                owner_ = true;
                return *this;
            }

            mat_ = std::move(other.mat_);
            owner_ = std::move(other.owner_);
            return *this;
        }

        void finalize();

        void fillComplete()
        {
            mat_->fillComplete();
        }

        void replaceGlobalValues (const GO globalRow, const LO numEnt, const SC vals[], const GO cols[])
        {
            mat_->replaceGlobalValues(globalRow, numEnt, vals, cols);
        }

        void replaceLocalValues (const LO localRow, const LO numEnt,  const SC vals[], const LO cols[] )
        {
            mat_->replaceLocalValues(localRow, numEnt, vals, cols);
        }

        rcp_comm_type communicator() const
        {
            return implementation().getMap()->getComm();
        }

        void set_owner(const bool owner)
        {
            owner_ = owner;
        }

        //API functions
        void crs_init(const rcp_comm_type &comm,
                      std::size_t rows_local,
                      std::size_t cols_local,
                      Tpetra::global_size_t rows_global,
                      Tpetra::global_size_t cols_global,
                      std::size_t nnz_x_row);

        void crs_identity(const rcp_comm_type &comm,
                      std::size_t rows_local,
                      std::size_t cols_local,
                      Tpetra::global_size_t rows_global,
                      Tpetra::global_size_t cols_global,
                      const Scalar factor = 1.);


        inline Range row_range() const
        {
            return  { implementation().getRowMap()->getMinGlobalIndex(), implementation().getRowMap()->getMaxGlobalIndex() + 1 };
        }

        inline Range col_range() const
        {
            if(implementation().getDomainMap().is_null()) {
                assert(!init_->domain_map.is_null());
                return  { init_->domain_map->getMinGlobalIndex(), init_->domain_map->getMaxGlobalIndex() + 1 };
            } else {
                return  { implementation().getDomainMap()->getMinGlobalIndex(), implementation().getDomainMap()->getMaxGlobalIndex() + 1 };
            }
        }

        inline Size size() const
        {
            if(is_null()) {
                return {0, 0};
            }

            if(implementation().isFillComplete()) {
                return { implementation().getGlobalNumRows(), implementation().getGlobalNumCols() };
            } else {
                assert(!implementation().getRowMap().is_null());
                assert(!implementation().getColMap().is_null());

                return { implementation().getRowMap()->getGlobalNumElements(), implementation().getColMap()->getGlobalNumElements() };
            }
        }

        inline Size local_size() const
        {
            assert(!implementation().getRowMap().is_null());

            if(implementation().getDomainMap().is_null()) {
                assert(!init_->domain_map.is_null());
                return { implementation().getRowMap()->getNodeNumElements(), init_->domain_map->getNodeNumElements() };
            } else {
                return { implementation().getRowMap()->getNodeNumElements(), implementation().getDomainMap()->getNodeNumElements() };
            }
        }

        inline void read_lock()
        {
            //TODO?
        }

        inline void read_unlock()
        {
            //TODO?
        }

        inline void write_lock()
        {
            //TODO?
            implementation().resumeFill();
        }

        inline void write_unlock()
        {
            this->finalize();
        }

        void describe(std::ostream &os) const
        {
            auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(os));
            implementation().describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
        }

        inline void describe() const
        {
            describe(std::cout);
        }

        void set(const GO &row, const GO &col, const Scalar &value);
        Scalar get(const GO &row, const GO &col) const;
        void add(const GO &row, const GO &col, const Scalar &value);

        template<typename Integer>
        void add_matrix(
         const std::vector<Integer> &rows,
         const std::vector<Integer> &cols,
         const std::vector<Scalar> &values
         )
        {
            //FIXME and find more efficient way
            const auto n_rows = rows.size();
            const auto n_cols = cols.size();
            assert(values.size() == n_rows*n_cols);

            for(std::size_t i = 0; i < n_rows; ++i) {
                const auto i_offset = i*n_rows;
                for(std::size_t j = 0; j < n_cols; ++j) {
                    add(rows[i], cols[j], values[i_offset + j]);
                }
            }
        }

       template<typename Integer>
        void set_matrix(
         const std::vector<Integer> &rows,
         const std::vector<Integer> &cols,
         const std::vector<Scalar> &values
         )
        {
            //FIXME and find more efficient way
           const auto n_rows = rows.size();
           const auto n_cols = cols.size();
           assert(values.size() == n_rows*n_cols);

           for(std::size_t i = 0; i < n_rows; ++i) {
               const auto i_offset = i*n_rows;
               for(std::size_t j = 0; j < n_cols; ++j) {
                   set(rows[i], cols[j], values[i_offset + j]);
               }
           }
        }

        void mult(const TpetraVector &vec, TpetraVector &result) const;
        void mult_t(const TpetraVector &vec, TpetraVector &result) const;

        void mult(const TpetraMatrix &right, TpetraMatrix &result) const;
        //result = tranpose(*this) * mat
        void mult_t(const TpetraMatrix &right, TpetraMatrix &result) const;

        //result op(*this) * op
        void mult(const bool transpose_this, const TpetraMatrix &right, const bool transpose_right, TpetraMatrix &result) const;
        void axpy(const Scalar alpha, const TpetraMatrix &x);
        void transpose(TpetraMatrix &mat) const;

        void get_diag(TpetraVector &d) const;
        void init_diag(const TpetraVector &d);

        inline void scale(const Scalar alpha)
        {
            //Why???
            write_lock();
            implementation().scale(alpha);
            write_unlock();
        }

        inline crs_mat_type &implementation()
        {
            assert(!mat_.is_null());
            return *mat_;
        }

        inline const crs_mat_type &implementation() const
        {
            assert(!mat_.is_null());
            return *mat_;
        }

        inline rcp_crs_mat_type &implementation_ptr()
        {
            assert(!mat_.is_null());
            return mat_;
        }

        inline const rcp_crs_mat_type &implementation_ptr() const
        {
            assert(!mat_.is_null());
            return mat_;
        }

        inline bool is_null() const
        {
            return mat_.is_null();
        }

        bool read(const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const std::string &path);
        bool write(const std::string &path) const;

        bool is_valid(const bool verbose = false) const;

        Scalar norm2() const;
        Scalar sum() const;


    private:
        rcp_crs_mat_type  mat_;
        bool              owner_;

        typedef struct {
            rcp_map_type domain_map;
            rcp_map_type range_map;
        } InitStructs;

        std::shared_ptr<InitStructs> init_;
    }; //TpetraMatrix
}

#endif //UTOPIA_TPETRAMATRIX_H

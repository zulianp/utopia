#ifndef UTOPIA_TPETRAMATRIX_H
#define UTOPIA_TPETRAMATRIX_H

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"
#include "utopia_Size.hpp"
#include "utopia_Logger.hpp"

#include "utopia_Tpetra_Vector.hpp"

#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_Map_decl.hpp>


#include <iostream>
#include <memory>

namespace utopia {

    class TpetraMatrix {
    public:

        typedef Tpetra::CrsMatrix<>                       crs_matrix_type;
        typedef Teuchos::RCP<crs_matrix_type>             rcp_crs_matrix_type;
        typedef Teuchos::RCP<const Teuchos::Comm<int> >   rcp_comm_type;
        typedef Tpetra::Map<>                             map_type;
        typedef Teuchos::RCP<map_type>                    rcp_map_type;
        typedef Tpetra::Vector<>::local_ordinal_type      local_ordinal_type;
        typedef Tpetra::Vector<>::global_ordinal_type     global_ordinal_type;
        typedef Tpetra::Vector<>::scalar_type             Scalar;
        typedef crs_matrix_type::node_type node_type;

        TpetraMatrix() : owner_(true) {}

        /////////////////////////////////////////////////////////////
        ~TpetraMatrix()
        {}

        //deep copy
        //     template <class Node2>
        // rcp_crs_matrix_type  clone (
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

        void set(const global_ordinal_type &row, const global_ordinal_type &col, const Scalar &value);
        Scalar get(const global_ordinal_type &row, const global_ordinal_type &col) const;
        void add(const global_ordinal_type &row, const global_ordinal_type &col, const Scalar &value);

        template<typename Integer>
        void add_matrix(
         const std::vector<Integer> &rows,
         const std::vector<Integer> &cols,
         const std::vector<Scalar> &values
         )
        {
            //FIXME and find more efficient way
            const auto n = values.size();
            assert(n == rows.size());
            assert(n == cols.size());

            for(std::size_t i = 0; i < n; ++i) {
                add(rows[i], cols[i], values[i]);
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
            const auto n = values.size();
            assert(n == rows.size());
            assert(n == cols.size());

            for(std::size_t i = 0; i < n; ++i) {
                set(rows[i], cols[i], values[i]);
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

        inline crs_matrix_type &implementation()
        {
            assert(!mat_.is_null());
            return *mat_;
        }

        inline const crs_matrix_type &implementation() const
        {
            assert(!mat_.is_null());
            return *mat_;
        }

        inline bool is_null() const
        {
            return mat_.is_null();
        }

        bool read(const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const std::string &path);
        bool write(const std::string &path) const;

        bool is_valid(const bool verbose = false) const;

        Scalar norm2() const;


    private:
        rcp_crs_matrix_type  mat_;
        bool                 owner_;

        typedef struct {
            rcp_map_type domain_map;
            rcp_map_type range_map;
        } InitStructs;

        std::shared_ptr<InitStructs> init_;
    }; //TpetraMatrix
}

#endif //UTOPIA_TPETRAMATRIX_H

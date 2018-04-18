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
        
        inline void finalize()
        {
            implementation().fillComplete();
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
        
        
        inline Range row_range() const
        {
            return  { implementation().getRowMap()->getMinGlobalIndex(), implementation().getRowMap()->getMaxGlobalIndex() + 1 };
        }
        
        inline Size size() const
        {
            return { implementation().getRowMap()->getGlobalNumElements(), implementation().getColMap()->getGlobalNumElements() };
        }
        
        inline Size local_size() const
        {
            return { implementation().getRowMap()->getNodeNumElements(), implementation().getColMap()->getNodeNumElements() };
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
        
        void describe(std::ostream &os = std::cout) const
        {
            auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(os));
            implementation().describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
        }
        
        void set(const global_ordinal_type &row, const global_ordinal_type &col, const Scalar &value);
        void add(const global_ordinal_type &row, const global_ordinal_type &col, const Scalar &value);

        void mult(const TpetraVector &vec, TpetraVector &result) const;
        void mult(const TpetraMatrix &right, TpetraMatrix &result) const;
        void axpy(const Scalar alpha, const TpetraMatrix &x);
       
        inline void scale(const Scalar alpha)
        {
            implementation().scale(alpha);
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
        
    private:
        rcp_crs_matrix_type  mat_;
        bool                 owner_;
    }; //TpetraMatrix
}

#endif //UTOPIA_TPETRAMATRIX_H

#ifndef UTOPIA_TPETRAMATRIX_H
#define UTOPIA_TPETRAMATRIX_H

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"
#include "utopia_Size.hpp"
#include "utopia_Matrix.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_Normed.hpp"
#include "utopia_Operator.hpp"

#include "utopia_Tpetra_Vector.hpp"
#include "utopia_trilinos_Communicator.hpp"

#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <iostream>
#include <memory>

namespace utopia {
    //template<class NodeType>
    class TpetraMatrix :
    public DistributedSparseMatrix<TpetraScalar, TpetraSizeType>,
    public SparseConstructible<TpetraScalar, TpetraSizeType>,
    public Normed<TpetraScalar>,
    public Transformable<TpetraScalar>,
    // Static polymorphic types
    public BLAS1Tensor<TpetraMatrix>,
    public BLAS2Matrix<TpetraMatrix, TpetraVector>,
    public BLAS3Matrix<TpetraMatrix>,
    public Comparable<TpetraMatrix>,
    public Operator<TpetraVector>,
    public Tensor<TpetraMatrix, 2>
    {
    public:

        /////////////////////////////////////////////////////////////
        // typedef definitions
        /////////////////////////////////////////////////////////////

        using Scalar   = utopia::TpetraScalar;
        using SizeType = utopia::TpetraSizeType;

        //types of Operators
        typedef Tpetra::Operator<>::scalar_type SC;
        typedef Tpetra::Operator<SC>::local_ordinal_type LO;
        typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;

        //types of Kokkos Parallel Nodes
        typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;
#ifdef KOKKOS_ENABLE_CUDA
        typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
        typedef cuda_node NT;
#elif defined KOKKOS_ENABLE_ROCM //Kokkos::Compat::KokkosROCmWrapperNode doesn't exist
        typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::ROCm> rocm_node;
        typedef rocm_node NT;
#elif defined KOKKOS_ENABLE_OPENMP
        typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
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

        using IndexSet    = Traits<TpetraMatrix>::IndexSet;
        using IndexArray  = Traits<TpetraMatrix>::IndexArray;
        using ScalarArray = Traits<TpetraMatrix>::ScalarArray;

        /////////////////////////////////////////////////////////////
        //Constructors
        /////////////////////////////////////////////////////////////

        //Default Constructor
        TpetraMatrix() : owner_(true) {}

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

        TpetraMatrix(const rcp_crs_mat_type &mat, const bool owner = false)
        : mat_(mat), owner_(owner)
        {}

        /////////////////////////////////////////////////////////////
        //Destructor
        /////////////////////////////////////////////////////////////

         ~TpetraMatrix()
         {}


        ////////////////////////////////////////////////////////////////////
        ///////////////////////// BOILERPLATE CODE FOR EDSL ////////////////
        ////////////////////////////////////////////////////////////////////

         using Super    = utopia::Tensor<TpetraMatrix, 2>;
         using Super::Super;

         template<class Expr>
         TpetraMatrix(const Expression<Expr> &expr)
        : owner_(true)
         {
            static_assert(!std::is_same<Expr, Tensor<TpetraVector, 1>>::value, "cannot assign a vector to a matrix");
             //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
             Super::construct_eval(expr.derived());
         }

         template<class Expr>
         inline TpetraMatrix &operator=(const Expression<Expr> &expr)
         {
             Super::assign_eval(expr.derived());
             return *this;
         }

         void assign(const TpetraMatrix &other) override
         {
             copy(other);
         }

         void assign(TpetraMatrix &&other) override
         {
            owner_ = std::move(other.owner_);
            comm_ = std::move(other.comm_);
            mat_ = std::move(other.mat_);
         }


        void copy(const TpetraMatrix &other) override
        {
            if(this == &other) return;

            if(other.is_null()) {
                mat_.reset();
                owner_ = true;
                return;
            }

            mat_ = other.mat_->clone(other.mat_->getNode());
            owner_ = true;
        }


        void select(
            const IndexSet &row_index,
            const IndexSet &col_index,
            TpetraMatrix &result) const override;

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR DistributedObject ////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        TrilinosCommunicator &comm() override
        {
            return comm_;
        }

        const TrilinosCommunicator &comm() const override
        {
            return comm_;
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR MatrixBase and DistributedMatrix /////////////
        ///////////////////////////////////////////////////////////////////////////

        void c_set(const SizeType &i, const SizeType &j, const Scalar &value) override;
        void c_add(const SizeType &i, const SizeType &j, const Scalar &value) override;

        SizeType rows() const override;
        SizeType cols() const override;
        SizeType local_rows() const override;
        SizeType local_cols() const override;

        inline Range row_range() const override
        {
            return  { implementation().getRowMap()->getMinGlobalIndex(), implementation().getRowMap()->getMaxGlobalIndex() + 1 };
        }

        inline Range col_range() const override
        {
            if(implementation().getDomainMap().is_null()) {
                assert(!init_->domain_map.is_null());
                return  { init_->domain_map->getMinGlobalIndex(), init_->domain_map->getMaxGlobalIndex() + 1 };
            } else {
                return  { implementation().getDomainMap()->getMinGlobalIndex(), implementation().getDomainMap()->getMaxGlobalIndex() + 1 };
            }
        }

        inline Size size() const override
        {
            if(is_null()) {
                return {0, 0};
            }

            if(implementation().isFillComplete()) {
                return { implementation().getGlobalNumRows(), implementation().getGlobalNumCols() };
            } else {
                assert(!implementation().getRowMap().is_null());

                if(implementation().getDomainMap().is_null()) {
                    assert(!init_->domain_map.is_null());
                    return { implementation().getRowMap()->getGlobalNumElements(), init_->domain_map->getGlobalNumElements() };
                } else {
                    return { implementation().getRowMap()->getGlobalNumElements(), implementation().getDomainMap()->getGlobalNumElements() };
                }
            }
        }

        inline bool empty() const override
        {
            return is_null();
        }

        inline Size local_size() const override
        {
            if(is_null()) {
                return {0, 0};
            }

            assert(!implementation().getRowMap().is_null());

            if(implementation().getDomainMap().is_null()) {
                assert(!init_->domain_map.is_null());
                return { implementation().getRowMap()->getNodeNumElements(), init_->domain_map->getNodeNumElements() };
            } else {
                return { implementation().getRowMap()->getNodeNumElements(), implementation().getDomainMap()->getNodeNumElements() };
            }
        }

        void describe(std::ostream &os) const
        {
            auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(os));
            implementation().describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
        }

        inline void describe() const override
        {
            describe(std::cout);
        }

        void clear() override;

        /////////////////////////////////////////////////////////////
        //OVERRIDES for SparseConstructible
        /////////////////////////////////////////////////////////////

        void identity(const Size &s, const Scalar &diag = 1.0) override;

        ///Specialize for sparse matrices
        void sparse(const Size &s, const SizeType &/*nnz*/) override;

        ///Specialize for sparse matrices
        void local_sparse(const Size &s, const SizeType &/*nnz*/) override;

        void local_identity(const Size &s, const Scalar &diag = 1.0) override;

        /////////////////////////////////////////////////////////////
        //OVERRIDES for Normed
        /////////////////////////////////////////////////////////////


        Scalar norm_infty() const override;
        Scalar norm1() const override;
        Scalar norm2() const override;

        /////////////////////////////////////////////////////////////
        //Overloading Operators
        /////////////////////////////////////////////////////////////

        TpetraMatrix &operator=(const TpetraMatrix &other)
        {
            copy(other);
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

        // void fillComplete()
        // {
        //     mat_->fillComplete();
        // }

        // void replaceGlobalValues (const GO globalRow, const LO numEnt, const SC vals[], const GO cols[])
        // {
        //     mat_->replaceGlobalValues(globalRow, numEnt, vals, cols);
        // }

        // void replaceLocalValues (const LO localRow, const LO numEnt,  const SC vals[], const LO cols[] )
        // {
        //     mat_->replaceLocalValues(localRow, numEnt, vals, cols);
        // }

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

        void crs_init(const rcp_comm_type &comm,
                      std::size_t rows_local,
                      std::size_t cols_local,
                      Tpetra::global_size_t rows_global,
                      Tpetra::global_size_t cols_global,
                      const Teuchos::ArrayRCP<size_t> &rowPtr,
                      const Teuchos::ArrayRCP<LO> &cols,
                      const Teuchos::ArrayRCP<Scalar> &values);

        void crs_identity(const rcp_comm_type &comm,
                      std::size_t rows_local,
                      std::size_t cols_local,
                      Tpetra::global_size_t rows_global,
                      Tpetra::global_size_t cols_global,
                      const Scalar factor);


        inline void read_lock() override
        {
            //TODO?
        }

        inline void read_unlock() override
        {
            //TODO?
        }

        inline void write_lock(WriteMode mode = utopia::AUTO) override
        {
            //TODO?
            implementation().resumeFill();
        }

        inline void write_unlock(WriteMode mode = utopia::AUTO) override
        {
            this->finalize();
        }

        void set(const SizeType &row, const SizeType &col, const Scalar &value) override;
        Scalar get(const SizeType &row, const SizeType &col) const;
        void add(const SizeType &row, const SizeType &col, const Scalar &value) override;

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

       template<typename Integer>
        void set_matrix(
         const std::vector<Integer> &rows,
         const std::vector<Integer> &cols,
         const Scalar value
         )
        {
            //FIXME and find more efficient way
           const auto n_rows = rows.size();
           const auto n_cols = cols.size();
           for(std::size_t i = 0; i < n_rows; ++i) {
               for(std::size_t j = 0; j < n_cols; ++j) {
                   set(rows[i], cols[j], value);
               }
           }
        }

        inline bool apply(const TpetraVector &in, TpetraVector &out) const override
        {
            if(empty()) return false;

            this->multiply(in, out);

            return true;
        }

        void multiply(const TpetraVector &vec, TpetraVector &result) const override;
        void transpose_multiply(const TpetraVector &vec, TpetraVector &result) const override;

        void multiply(const TpetraMatrix &right, TpetraMatrix &result) const override;
        void multiply(const Scalar &alpha, const TpetraMatrix &B, TpetraMatrix &C) const override;
        //result = tranpose(*this) * mat
        void transpose_multiply(const TpetraMatrix &right, TpetraMatrix &result) const override;

        void multiply(
            const bool transpose_A,
            const bool transpose_B,
            const TpetraMatrix &B,
            TpetraMatrix &C) const override;
        
        void axpy(const Scalar &alpha, const TpetraMatrix &x) override;
        void transpose(TpetraMatrix &mat) const override;

        void build_diag(TpetraVector &d) const;
        void diag(const TpetraVector &d);
        void diag(const TpetraMatrix &d);



        inline rcp_crs_mat_type &raw_type()
        {
           return implementation_ptr();
        }

        inline const rcp_crs_mat_type &raw_type() const
        {
           return implementation_ptr();
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

        inline bool read(const std::string &path)
        {
            return read(comm().get(), path);
        }

        bool read(const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const std::string &path);
        bool write(const std::string &path) const;

        bool is_valid(const bool verbose = false) const;

        Scalar sum() const;
        //FIXME
        inline Scalar reduce(const Plus &) const //override
        {
            return sum();
        }


        void set_domain_and_range(
            const rcp_map_type &domain_map,
            const rcp_map_type &range_map
            )
        {
            init_ = std::make_shared<InitStructs>();
            init_->domain_map = domain_map;
            init_->range_map = range_map;
        }

        void build_from_structure(const TpetraMatrix &rhs);


        void set_zero_rows(const IndexSet &index, const Scalar &diag);

        ////////////////////////////////////////////////////////////////////////
        //////////////////////////// Transormable //////////////////////////////
        ////////////////////////////////////////////////////////////////////////


        void transform(const Sqrt &) override;
        void transform(const Pow2 &) override;
        void transform(const Log &) override;
        void transform(const Exp &) override;
        void transform(const Cos &) override;
        void transform(const Sin &) override;
        void transform(const Abs &) override;
        void transform(const Minus &) override;

        void transform(const Pow &p) override;
        void transform(const Reciprocal<Scalar> &f) override;

        ////////////////////////////////////////////////////////////////////////
        //////////////////////////// Blas1Tensor //////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        void swap(TpetraMatrix &x) override;
        void scale(const Scalar &alpha) override;
        Scalar dot(const TpetraMatrix &other) const override;


        ////////////////////////////////////////////////////////////////////////
        //////////////////////////// Blas2Matrix //////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        void gemv(const bool transpose, const Scalar &alpha, const TpetraVector &x, const Scalar &beta, TpetraVector &y) const override;

        ////////////////////////////////////////////////////////////////////////
        //////////////////////////// Comparable //////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        bool equals(const TpetraMatrix &other, const Scalar &tol = 0.0) const override;



        inline std::string get_class() const override
        {
            return "TpetraMatrix";
        }
        
        

    private:
        TrilinosCommunicator comm_;
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

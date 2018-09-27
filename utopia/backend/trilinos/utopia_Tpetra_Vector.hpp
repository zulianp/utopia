
#ifndef UTOPIA_TPETRA_VECTOR_HPP
#define UTOPIA_TPETRA_VECTOR_HPP

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"
#include "utopia_Size.hpp"
#include "utopia_Writable.hpp"

#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Vector_decl.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>


#include <memory>

namespace utopia {

    class TpetraVector {
    public: 
    
    
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

    typedef Tpetra::Operator<SC, LO, GO, NT> OP;

    typedef Tpetra::Map<LO, GO, NT>                   map_type;
    typedef Tpetra::Vector<SC, LO, GO, NT>            vector_type;
    typedef Teuchos::RCP<vector_type>                 rcpvector_type;
    typedef Teuchos::RCP<const Teuchos::Comm<int> >   rcp_comm_type;
    typedef Teuchos::RCP<const map_type>              rcp_map_type;

//        typedef Tpetra::Vector<>::mag_type                magnitude_type;
    typedef vector_type::scalar_type                  Scalar;

        TpetraVector()
        {
            int indexBase = 0;
            auto comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
            auto conigMap = Teuchos::rcp (new map_type (comm->getSize (), indexBase, comm));
            vec_.reset(new vector_type (conigMap, false));
        }

        ~TpetraVector()
        { }

        TpetraVector(const TpetraVector &other);


        TpetraVector(TpetraVector &&other)
        : vec_(std::move(other.vec_))
        { }

        rcp_comm_type communicator()
        {
            return implementation().getMap()->getComm();
        }

        const rcp_comm_type communicator() const
        {
            return implementation().getMap()->getComm();
        }

        TpetraVector &operator=(const TpetraVector &other);

        TpetraVector &operator=(TpetraVector &&other)
        {
            if(this == &other) return *this;
            vec_ = std::move(other.vec_);
            ghosted_vec_ = std::move(other.ghosted_vec_);
            return *this;
        }

        void copy(const TpetraVector &other);


        //////////////////////////////////////////
        //API functions
        //////////////////////////////////////////
        inline void values(const rcp_comm_type &comm, std::size_t n_local, Tpetra::global_size_t n_global, Scalar value)
        {
            rcp_map_type map;

            if(n_local == INVALID_INDEX) {
                map = Teuchos::rcp(new map_type(n_global, 0, comm));
            } else {
                map = Teuchos::rcp(new map_type(n_global, n_local, 0, comm));
            }

            vec_.reset(new vector_type(map));
            implementation().putScalar(value);
        }

        inline void init(const rcp_map_type &map)
        {
            vec_.reset(new vector_type(map));
        }



        void ghosted(const rcp_comm_type &comm, 
                     const TpetraVector::GO &local_size,
                     const TpetraVector::GO &global_size,
                     const std::vector<GO> &ghost_index
        );

        inline void axpy(const Scalar &alpha, const TpetraVector &x)
        {
            implementation().update(alpha, *x.vec_, 1.);
        }

        void describe(std::ostream &os) const
        {
            auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(os));
            implementation().describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
        }


        void describe() const
        {
            describe(std::cout);
        }

        inline Scalar get(const GO i) const
        {
            assert(!read_only_data_.is_null() && "Use Read<Vector> w(v); to enable reading from this vector v!");
            return read_only_data_[local_index(i)];

        }

        inline LO local_index(const GO i) const
        {
            if(has_ghosts()) {
               return ghosted_vec_->getMap()->getLocalElement(i);
            } else {
                //i - implementation().getMap()->getMinGlobalIndex()
                return implementation().getMap()->getLocalElement(i);
            }
        }

        inline void set(const GO i, const Scalar value)
        {
            if(!write_data_.is_null()) {
                write_data_[i - implementation().getMap()->getMinGlobalIndex()] = value;
            } else {
                implementation().replaceGlobalValue(i, value);
            }
        }

        inline void add(const GO i, const Scalar value)
        {
            if(!ghosted_vec_.is_null()) {
                ghosted_vec_->sumIntoGlobalValue(i, value);
            } else {
                implementation().sumIntoGlobalValue(i, value);
            }
        }

        inline void set(const Scalar value)
        {
            implementation().putScalar(value);
        }

        template<typename Integer>
        void get(
            const std::vector<Integer> &index,
            std::vector<Scalar> &values) const
        {
            // m_utopia_warning_once(" > get does not work in parallel if it asks for ghost entries");
            //FIXME does not work in parallel
            auto data  = get_read_only_data();
            // auto offset = implementation().getMap()->getMinGlobalIndex();

            auto n = index.size();
            values.resize(n);

            for(std::size_t i = 0; i < n; ++i) {
                auto local_index = this->local_index(index[i]);// - offset;
                assert(local_index < data.size());
                values[i] = data[local_index];
            }
        }

        void set_vector(
            const std::vector<GO> &indices,
            const std::vector<Scalar> &values);

        void add_vector(
            const std::vector<GO> &indices,
            const std::vector<Scalar> &values);

        template<typename Integer>
        void select(
            const std::vector<Integer> &index,
            TpetraVector &out
            ) const
        {
            //FIXME does not work in parallel
            auto data   = implementation().getData();
            auto offset = implementation().getMap()->getMinGlobalIndex();

            auto map = Teuchos::rcp(
                new map_type(
                    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                    index.size(),
                    0,
                    communicator())
            );


            out.init(map);

            // if(communicator()->getSize() == 1) {

                auto n = data.size();
                auto out_data = out.implementation().getDataNonConst();

                for(std::size_t i = 0; i < n; ++i) {
                    auto local_index = index[i] - offset;
                    assert(local_index < n);
                    out_data[i] = data[local_index];
                }
            // } else {
            //     /////////////////////////////////////////////////

            //     std::vector<GO> tpetra_index;
            //     tpetra_index.reserve(index.size());

            //     for(auto i : index) {
            //         tpetra_index.push_back(i);
            //     }

            //     const Teuchos::ArrayView<const GO>
            //        index_view(tpetra_index);

            //      auto import_map = Teuchos::rcp(new map_type(global_size, index_view, 0, comm));

            //     Tpetra::Import<
            //         LO,
            //         GO,
            //         vector_type::node_type> importer(map, import_map);

            //     implementation().doImport(out.implementation(), importer, Tpetra::INSERT);
            // }
        }

        inline Teuchos::ArrayRCP<const Scalar> get_read_only_data() const
        {
            if(!ghosted_vec_.is_null()) {
                return ghosted_vec_->getData();
            } else {
                return implementation().getData();
            }
        }


        inline void read_lock()
        {
            read_only_data_ = get_read_only_data();

        }

        inline void read_unlock()
        {
            read_only_data_ = Teuchos::ArrayRCP<const Scalar>();
        }

        inline void write_lock(WriteMode mode)
        {
            //TODO?
        }

        void write_unlock(WriteMode mode);

        inline void read_and_write_lock()
        {
            write_data_ = implementation().getDataNonConst();
            read_only_data_ = write_data_;
        }

        inline void read_and_write_unlock()
        {
            write_data_ = Teuchos::ArrayRCP<Scalar>();
            read_only_data_ = Teuchos::ArrayRCP<const Scalar>();
        }

        inline Range range() const
        {
            return { implementation().getMap()->getMinGlobalIndex(), implementation().getMap()->getMaxGlobalIndex() + 1 };
        }

        inline Size size() const
        {
            if(is_null()) {
                return {INVALID_INDEX};
            }

            return { implementation().getMap()->getGlobalNumElements() };
        }

        inline Size local_size() const
        {
            if(is_null()) {
                return {INVALID_INDEX};
            }

            return { implementation().getMap()->getNodeNumElements() };
        }

        inline Scalar norm2() const {
           return implementation().norm2();
        }

        inline Scalar norm1() const {
           return implementation().norm1();
        }

        inline Scalar norm_infty() const {
          return implementation().normInf();
        }

        Scalar sum() const;
        Scalar min() const;
        Scalar max() const;

        bool is_nan_or_inf() const;

        inline void scale(const Scalar alpha)
        {
            implementation().scale(alpha);
        }

        inline Scalar dot(const TpetraVector &other) const
        {
            return implementation().dot(other.implementation());
        }


        inline void e_mul(const TpetraVector &right, TpetraVector &result) const
        {
            //https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1MultiVector.html#a95fae4b1f2891d8438b7fb692a85b3bd
            result.values(this->communicator(), local_size().get(0), size().get(0), 0.);
          
            result.implementation().elementWiseMultiply(
                1.,
                this->implementation(),
                right.implementation(),
                0.
            );
        }

        template<typename Op>
        inline void apply(const Op op)
        {
            read_and_write_lock();

            assert(write_data_.size() > 0);

            for(auto i = 0; i < write_data_.size(); ++i) {
                write_data_[i] = op.apply(write_data_[i]);
            }

            read_and_write_unlock();
        }

        void reciprocal(TpetraVector &result) const
        {
            if(result.empty() || result.size() != this->size())
            {
                result.init(this->implementation().getMap());
            }

            result.implementation().reciprocal(this->implementation());
        }

        template<typename Op>
        inline void apply_binary(const Op op, const TpetraVector &rhs, TpetraVector &result) const
        {
            assert(!empty());
            assert(!rhs.empty());
            assert(rhs.size() == size());
            assert(rhs.local_size() == local_size());

            if(result.empty() || result.size() != rhs.size())
            {
                result.init(rhs.implementation().getMap());
            } 

            auto a_lhs = this->implementation().getData();
            auto a_rhs = rhs.implementation().getData();
            auto a_res = result.implementation().getDataNonConst();

            assert(a_res.size() == a_lhs.size());
            assert(a_res.size() == a_rhs.size());

            for(auto i = 0; i < a_lhs.size(); ++i) {
               a_res[i] = op.apply(a_lhs[i], a_rhs[i]);
            }
        }

        inline vector_type &implementation()
        {
            assert(!vec_.is_null());
            return *vec_;
        }

        inline const vector_type &implementation() const
        {
            assert(!vec_.is_null());
            return *vec_;
        }

        inline rcpvector_type &implementation_ptr()
        {
            assert(!vec_.is_null());
            return vec_;
        }

        inline const rcpvector_type &implementation_ptr() const
        {
            assert(!vec_.is_null());
            return vec_;
        }

        inline bool is_null() const
        {
            return vec_.is_null();
        }

        void replaceGlobalValue (const GO globalRow, const SC &value )
        {
            std::cout << " sono qui " << std::endl;
            vec_->replaceGlobalValue (globalRow, value);

        }

        void replaceLocalValue (const LO localRow, const SC &value )
        {
            std::cout << localRow << " localRow " << std::endl;
            std::cout << value << " value " << std::endl;
            vec_->replaceLocalValue(localRow, value);
        }

        bool read(const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const std::string &path);
        bool write(const std::string &path) const;

        inline bool empty() const
        {
            return vec_.is_null();
        }

        void update_ghosts();
        void export_ghosts_add();

        bool has_ghosts() const
        {
            return !ghosted_vec_.is_null();
        }
        
    private:
        rcpvector_type vec_;
        rcpvector_type ghosted_vec_;
        Teuchos::ArrayRCP<const Scalar> read_only_data_;
        Teuchos::ArrayRCP<Scalar> write_data_;
    };
}

#endif //UTOPIA_TPETRA_VECTOR_HPP


#ifndef UTOPIA_TPETRA_VECTOR_HPP
#define UTOPIA_TPETRA_VECTOR_HPP

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"
#include "utopia_Size.hpp"

#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Vector_decl.hpp>

#include <memory>

namespace utopia {

    class TpetraVector {
    public:
        typedef Tpetra::Map<>                             map_type;
        typedef Tpetra::Vector<>                          vector_type;
        typedef Teuchos::RCP<vector_type>                 rcpvector_type;
        typedef Teuchos::RCP<const Teuchos::Comm<int> >   rcp_comm_type;
        typedef Teuchos::RCP<const map_type>              rcp_map_type;
        typedef Tpetra::Vector<>::scalar_type             scalar_type;
        typedef Tpetra::Vector<>::local_ordinal_type      local_ordinal_type;
        typedef vector_type::global_ordinal_type          global_ordinal_type;
        typedef Tpetra::Vector<>::mag_type                magnitude_type;
        typedef vector_type::scalar_type Scalar;

        TpetraVector()
        {
            // FIXME global size to size of the comm and index base to zero
            // auto comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
            // auto conigMap = Teuchos::rcp (new Tpetra::Map<> (comm->getSize (), 0, comm));
            // vec_.reset(new vector_type (conigMap));
        }

        ~TpetraVector()
        { }

        TpetraVector(const TpetraVector &other)
        : vec_(Teuchos::rcp(new vector_type(*other.vec_, Teuchos::Copy)))
        { }


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

        TpetraVector &operator=(const TpetraVector &other)
        {
            if(this == &other) return *this;

            if(other.is_null()) {
                vec_.reset();
                return *this;
            }

            vec_ = Teuchos::rcp(new vector_type(*other.vec_, Teuchos::Copy));
            return *this;
        }

        TpetraVector &operator=(TpetraVector &&other)
        {
            if(this == &other) return *this;
            vec_ = std::move(other.vec_);
            return *this;
        }


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

            vec_ = Teuchos::rcp(new vector_type(map));
            implementation().putScalar(value);
        }

        inline void init(const rcp_map_type &map)
        {
            vec_ = Teuchos::rcp(new vector_type(map));
        }

        inline void axpy(const Scalar &alpha, const TpetraVector &x)
        {
            implementation().update(alpha, *x.vec_, 1.);
        }

        inline void describe(std::ostream &os) const
        {
            auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(os));
            implementation().describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
        }


        void describe() const
        {
            describe(std::cout);
        }

        inline Scalar get(const global_ordinal_type i) const
        {
            assert(!read_only_data_.is_null() && "Use Read<Vector> w(v); to enable reading from this vector v!");
            return read_only_data_[i - implementation().getMap()->getMinGlobalIndex()];
        }

        inline void set(const global_ordinal_type i, const Scalar value)
        {
            if(!write_data_.is_null()) {
                write_data_[i - implementation().getMap()->getMinGlobalIndex()] = value;
            } else {
                implementation().replaceGlobalValue(i, value);
            }
        }

        inline void add(const global_ordinal_type i, const Scalar value)
        {
            implementation().sumIntoGlobalValue(i, value);
        }

        inline void set(const Scalar value)
        {
            implementation().putScalar(value);
        }

        inline void read_lock()
        {
            read_only_data_ = implementation().getData();
        }

        inline void read_unlock()
        {
            read_only_data_ = Teuchos::ArrayRCP<const Scalar>();
        }

        inline void write_lock()
        {
            //TODO?
        }

        inline void write_unlock()
        {
            //TODO?
        }

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

        inline Scalar sum() const {
          // return implementation().sum();
            assert(false && "implement me");
            return -1.;
        }

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
                result.implementation(),
                right.implementation(),
                0.
            );
        }

        inline vector_type &implementation()
        {
            return *vec_;
        }

        inline const vector_type &implementation() const
        {
            return *vec_;
        }

        inline bool is_null() const
        {
            return vec_.is_null();
        }

        bool read(const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const std::string &path);
        bool write(const std::string &path) const;

    private:
        rcpvector_type vec_;
        Teuchos::ArrayRCP<const Scalar> read_only_data_;
        Teuchos::ArrayRCP<Scalar> write_data_;
    };
}

#endif //UTOPIA_TPETRA_VECTOR_HPP

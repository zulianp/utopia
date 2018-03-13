
#ifndef UTOPIA_TPETRAVECTOR_H
#define UTOPIA_TPETRAVECTOR_H

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"

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
        : vec_(Teuchos::rcp(new vector_type(*other.vec_)))
        { }

        rcp_comm_type communicator()
        {
            return vec_->getMap()->getComm();
        }

        const rcp_comm_type communicator() const
        {
            return vec_->getMap()->getComm();
        }

        TpetraVector &operator=(const TpetraVector &other)
        {
            if(this == &other) return *this;
            vec_ = Teuchos::rcp(new vector_type(*other.vec_));
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
            auto map = Teuchos::rcp(new map_type(n_global, n_local, 0, comm));
            vec_ = Teuchos::rcp(new vector_type(map));
            vec_->putScalar(value);
        }

        inline void axpy(const Scalar &alpha, const TpetraVector &x)
        {
            vec_->update(alpha, *x.vec_, 1.);
        }

        inline void describe(std::ostream &os = std::cout) const
        {
            auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(os));
            vec_->describe(*out,  Teuchos::EVerbosityLevel::VERB_EXTREME);
        }

    private:
        rcpvector_type vec_;
    };
}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H

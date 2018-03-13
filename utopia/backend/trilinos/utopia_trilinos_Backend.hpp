#ifndef UTOPIA_TRILINOS_BACKEND_HPP
#define UTOPIA_TRILINOS_BACKEND_HPP

#include "utopia_trilinos_Traits.hpp"
#include "utopia_Core.hpp"
#include "utopia_Factory.hpp"
//#include "utopia_Base.hpp"
#include "utopia_ScalarBackend.hpp"

#include <utility>

namespace utopia {
    class TrilinosBackend : public ScalarBackend<double> {
    public:
        typedef double Scalar;
        typedef TpetraVector Vector;
        typedef TpetraMatrix Matrix;
        typedef TpetraSparseMatrix SparseMatrix;

        static Range range(const TpetraVector &v)
        {
            return v.range();
        }

        static Range row_range(const TpetraMatrix &m);
        static Range col_range(const TpetraMatrix &m);

        static void size(const TpetraVector &v, Size &size)
        {
            size = v.size();
        }

        static void local_size(const TpetraVector &v, Size &size)
        {
            size = v.local_size();
        }

        // static void size(const TpetraMatrix &m, Size &size); 
        // static void local_size(const TpetraMatrix &m, Size &size);


        //[io]
        // read matrix
        static bool read(const std::string &path, TpetraMatrix &m);
        // write matrix
        static bool write(const std::string &path, const TpetraMatrix &m);

        // read vector
        static bool read(const std::string &path, TpetraVector &v);

        // write vector
        static bool write(const std::string &path, const TpetraVector &v);

        //[builders]
        inline static void build(TpetraVector &v, const Size &size, const LocalValues<Scalar> &values)
        {
            v.values(default_communicator(), size.get(0), Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), values.value());
        }

        inline static void set(TpetraVector &v, const int index, Scalar value)
        {
            v.set(index, value);
        }

        inline static void add(TpetraVector &v, const int index, Scalar value)
        {
            v.add(index, value);
        }

        //[host/device locks]
        template<class Tensor>
        static void read_lock(const Tensor &) {}

        template<class Tensor>
        static void read_unlock(const Tensor &) {}

        static void write_lock(TpetraVector &vec)
        {
            vec.write_lock();
        }

        static void write_unlock(TpetraVector &vec)
        {
            vec.write_unlock();
        }

        static void write_lock(const TpetraMatrix &mat);
        // {
        //     mat.write_lock();
        // }

        static void write_unlock(const TpetraMatrix &mat);
        // {
        //     mat.write_unlock();
        // }

        // reductions
        // static Scalar norm2(const TpetraMatrix &m);
        inline static Scalar norm2(const TpetraVector &v)
        {
            return v.norm2();
        }

        inline static Scalar norm1(const TpetraVector &v)
        {
            return v.norm1();
        }

        inline static Scalar norm_infty(const TpetraVector &v)
        {
            return v.norm_infty();
        }

        //blas 1

        inline static void axpy(TpetraVector &y, const Scalar alpha, const TpetraVector &x)
        {
            y.axpy(alpha, x);
        }

    private:

        inline static auto default_communicator() -> decltype( Tpetra::DefaultPlatform::getDefaultPlatform().getComm() )
        {
            return Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
        }
    };

    template<>
    class Backend< double , TRILINOS > : public TrilinosBackend {
    public:
        inline static Backend &Instance()
        {
            static Backend instance;
            return instance;
        }

        BackendInfo &info()
        {
            return info_;
        }

    private:
        BackendInfo info_;

        Backend()
        {
            info_.set_name("trilinos");
        }
    };
}

#endif //UTOPIA_TRILINOS_BACKEND_HPP

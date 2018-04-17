#ifndef UTOPIA_TRILINOS_BACKEND_HPP
#define UTOPIA_TRILINOS_BACKEND_HPP

#include "utopia_trilinos_Traits.hpp"
#include "utopia_Core.hpp"
#include "utopia_Factory.hpp"
//#include "utopia_Base.hpp"
#include "utopia_ScalarBackend.hpp"

#include <utility>

namespace utopia {
    class TrilinosBackend : public ScalarBackend<TpetraVector::Scalar> {
    public:
        typedef TpetraVector::Scalar Scalar;
        typedef TpetraVector Vector;
        typedef TpetraMatrix Matrix;
        typedef TpetraSparseMatrix SparseMatrix;

        using ScalarBackend<Scalar>::apply_binary;
        using ScalarBackend<Scalar>::axpy;

        template<class LorRValueMatrix>
        static void assign(TpetraMatrix &left, LorRValueMatrix &&right)
        {
            left = std::forward<LorRValueMatrix>(right);
        }
        
        template<class LorRValueVector>
        static void assign(TpetraVector &left, LorRValueVector &&right)
        {
            left = std::forward<LorRValueVector>(right);
        }

        static Range range(const TpetraVector &v)
        {
            return v.range();
        }

        static Range row_range(const TpetraMatrix &m)
        {
            return m.row_range();
        }

        // static Range col_range(const TpetraMatrix &m)
        // {
        //     return m.col_range();
        // }

        static void size(const TpetraVector &v, Size &size)
        {
            size = v.size();
        }

        static void local_size(const TpetraVector &v, Size &size)
        {
            size = v.local_size();
        }

        static void size(const TpetraMatrix &m, Size &size)
        {
            size = m.size();
        }

        static void local_size(const TpetraMatrix &m, Size &size)
        {
            size = m.local_size();
        }

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

        inline static void build(TpetraVector &v, const Size &size, const Values<Scalar> &values)
        {
            v.values(default_communicator(), INVALID_INDEX, size.get(0), values.value());
        }

        inline static void build(TpetraVector &v, const Size &size, const LocalZeros &)
        {
            m_utopia_status_once("> Build zeros is using build values");
            v.values(default_communicator(), size.get(0), Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), 0.);
        }

        inline static void build(TpetraVector &v, const Size &size, const Zeros &)
        {
            m_utopia_status_once("> Build zeros is using build values");
            v.values(default_communicator(), INVALID_INDEX, size.get(0), 0.);
        }

        static void build(TpetraSparseMatrix &m, const Size &size, const LocalNNZ<std::size_t> &nnz)
        {
            m.crs_init(default_communicator(),
              size.get(0),
              size.get(1),
              Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
              Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
              nnz.nnz());
        }

        static void build(TpetraSparseMatrix &m, const Size &size, const NNZ<std::size_t> &nnz)
        {
            m.crs_init(default_communicator(),
              INVALID_INDEX,
              INVALID_INDEX,
              size.get(0),
              size.get(1),
              nnz.nnz());
        }

        inline static void set(TpetraVector &v, const int index, Scalar value)
        {
            v.set(index, value);
        }

        inline static void add(TpetraVector &v, const int index, Scalar value)
        {
            v.add(index, value);
        }

        inline static void set(TpetraMatrix &m, const TpetraMatrix::global_ordinal_type &row, const TpetraMatrix::global_ordinal_type &col, const Scalar &value)
        {
            m.set(row, col, value);
        }

        inline static void add(TpetraMatrix &m, const TpetraMatrix::global_ordinal_type &row, const TpetraMatrix::global_ordinal_type &col, const Scalar &value)
        {
            m.add(row, col, value);
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

        static void write_lock(TpetraMatrix &mat)
        {
            mat.write_lock();
        }

        static void write_unlock(TpetraMatrix &mat)
        {
            mat.write_unlock();
        }

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

        inline static void scale(TpetraVector &x, const Scalar alpha)
        {
            x.scale(alpha);
        }

        inline static void scale(TpetraMatrix &x, const Scalar alpha)
        {
            x.scale(alpha);
        }

        inline static Scalar dot(const TpetraVector &x, const TpetraVector &y)
        {
            return x.dot(y);
        }


        //blas 2
        static void multiply(
            TpetraVector &result,
            bool transpose_left,
            const TpetraMatrix &left,
            bool transpose_right,
            const Vector &right)
        {
            assert(!transpose_right);
            assert(!transpose_left);
            //TODO implement transpoe left

            left.mult(right, result);
        }

        inline static void apply_binary(TpetraVector &result, const TpetraMatrix &left, const Multiplies &, const TpetraVector &right)
        {
            left.mult(right, result);
        }

        inline static void apply_binary(TpetraMatrix &result, const TpetraMatrix &left, const Multiplies &, const TpetraMatrix &right)
        {
            left.mult(right, result);
        }

        inline static void apply_binary(TpetraVector &result, const TpetraVector &left, const Plus &, const TpetraVector &right)
        {
            result = left;
            result.axpy(1., right);
        }

        inline static void apply_binary(TpetraVector &result, const TpetraVector &left, const Minus &, const TpetraVector &right)
        {
            result = left;
            result.axpy(-1., right);
        }

        inline static void apply_binary(TpetraMatrix &result, const TpetraMatrix &left, const Plus &, const TpetraMatrix &right)
        {
            result = left;
            result.axpy(1., right);
        }

        inline static void apply_binary(TpetraMatrix &result, const TpetraMatrix &left, const Minus &, const TpetraMatrix &right)
        {
            result = left;
            result.axpy(-1., right);
        }

        static void apply_binary(TpetraMatrix &result, const Scalar factor, const Multiplies &, const TpetraMatrix &mat)
        {
            result = mat;
            result.scale(factor);
        }

        static void apply_binary(TpetraVector &result, const Scalar factor, const Multiplies &, const TpetraVector &vec)
        {
            result = vec;
            result.scale(factor);
        }

        inline static bool compare(const TpetraVector &left, const TpetraVector &right, const ApproxEqual &comp) {
            TpetraVector diff;
            apply_binary(diff, left, Minus(), right);
            return norm_infty(diff) <= comp.tol();
        }
        

        // monitoring functions for iterative solvers (Cyrill)
        // UTOPIA_DEPRECATED_MSG("Remove me")
        template<class Tensor>
        static void monitor(const long &, Tensor &)
        {
            m_utopia_status_once("Remove the monitor interface plz");
        }

    private:

        inline static auto default_communicator() -> decltype( Tpetra::DefaultPlatform::getDefaultPlatform().getComm() )
        {
            return Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
        }
    };

    template<>
    class Backend<TpetraVector::Scalar, TRILINOS> : public TrilinosBackend {
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

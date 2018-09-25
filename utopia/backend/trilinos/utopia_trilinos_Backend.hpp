#ifndef UTOPIA_TRILINOS_BACKEND_HPP
#define UTOPIA_TRILINOS_BACKEND_HPP

#include "utopia_trilinos_Traits.hpp"
#include "utopia_Core.hpp"
#include "utopia_Factory.hpp"
//#include "utopia_Base.hpp"
#include "utopia_ScalarBackend.hpp"

#include <utility>
#include <cmath>
#include <TpetraExt_MatrixMatrix_def.hpp>

//TODO find the configuration for having this
// #include "TpetraExt_TripleMatrixMultiply_def.hpp"


//useful links:
//https://trilinos.org/docs/dev/packages/tpetra/doc/html/namespaceTpetra_1_1TripleMatrixMultiply.html
//see MultiplyRAP

namespace utopia {
    class TrilinosBackend : public ScalarBackend<TpetraVector::Scalar> {
    public:
        typedef TpetraVector::Scalar Scalar;
        typedef TpetraVector Vector;
        typedef TpetraMatrix Matrix;
        typedef TpetraSparseMatrix SparseMatrix;

        using ScalarBackend<Scalar>::apply_binary;
        using ScalarBackend<Scalar>::axpy;
        using ScalarBackend<Scalar>::assign;

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


        template<class T>
        auto raw_type(T &t) -> decltype(t.implementation_ptr())
        {
            return t.implementation_ptr();
        }

        template<class T>
        auto raw_type(const T &t) -> decltype(t.implementation_ptr())
        {
            return t.implementation_ptr();
        }

        static Range range(const TpetraVector &v)
        {
            return v.range();
        }

        static Range row_range(const TpetraMatrix &m)
        {
            return m.row_range();
        }

        static Range col_range(const TpetraMatrix &m)
        {
            return m.col_range();
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
        static bool read(const std::string &path, TpetraMatrix &m)
        {
            return m.read(default_communicator(), path);
        }
        // write matrix
        static bool write(const std::string &path, const TpetraMatrix &m)
        {
            return m.write(path);
        }

        // read vector
        static bool read(const std::string &path, TpetraVector &v)
        {
            return v.read(default_communicator(), path);
        }

        // write vector
        static bool write(const std::string &path, const TpetraVector &v)
        {
            return v.write(path);
        }

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
            // m_utopia_status_once("> Build zeros is using build values");
            v.values(default_communicator(), size.get(0), Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), 0.);
        }

        inline static void build(TpetraVector &v, const Size &size, const Zeros &)
        {
            // m_utopia_status_once("> Build zeros is using build values");
            v.values(default_communicator(), INVALID_INDEX, size.get(0), 0.);
        }

        static void build(TpetraMatrix &m, const Size &size, const LocalZeros &)
        {
            assert(false && "implement me");
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

        static void build(TpetraSparseMatrix &m, const Size &size, const LocalIdentity &)
        {
            m.crs_identity(default_communicator(),
              size.get(0),
              size.get(1),
              Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
              Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()
            );
        }

        static void build(TpetraSparseMatrix &m, const Size &size, const Identity &)
        {
            m.crs_identity(default_communicator(),
              INVALID_INDEX,
              INVALID_INDEX,
              size.get(0),
              size.get(1)
            );
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

        inline static void build(TpetraSparseMatrix &, const Size &, const Zeros &)
        {
            m_utopia_error("> Build zeros is using build values");
            assert(false);
        }

        template<class Integer>
        void build_ghosts(
            const TpetraVector::global_ordinal_type &local_size,
            const TpetraVector::global_ordinal_type &global_size,
            const std::vector<Integer> &index,
            TpetraVector &vec)
        {
            assert(default_communicator()->getSize() == 1 && "implement me: does not work in parallel yet");
            vec.values(default_communicator(), local_size, global_size, 0.);
        }

        inline static void build_from_structure(TpetraSparseMatrix &lhs, const TpetraSparseMatrix &rhs)
        {
            //IMPLEMENT ME
            lhs = rhs;
        }

        template<class Expr>
        static void build_blocks(TpetraMatrix &left, const Blocks<Expr> &blocks)
        {
            assert(false && "IMPLEMENT ME");
        }


        template<class Expr>
        static void build_blocks(TpetraVector &left, const Blocks<Expr> &blocks)
        {
            assert(false && "IMPLEMENT ME");
        }

        static Scalar get(const TpetraVector &v, const TpetraVector::global_ordinal_type &index)
        {
            return v.get(index);
        }

        template<typename Integer>
       	static void get(const TpetraVector &v, const std::vector<Integer> &index, std::vector<Scalar> &values)
       	{
       		v.get(index, values);
       	}

        inline static void set(TpetraVector &v, const TpetraVector::global_ordinal_type &index, Scalar value)
        {
            v.set(index, value);
        }

        inline static void add(TpetraVector &v, const TpetraVector::global_ordinal_type &index, Scalar value)
        {
            v.add(index, value);
        }

        static Scalar get(const TpetraMatrix &m, const TpetraMatrix::global_ordinal_type &row, const TpetraMatrix::global_ordinal_type &col)
        {
            return m.get(row, col);
        }

        inline static void set(TpetraMatrix &m, const TpetraMatrix::global_ordinal_type &row, const TpetraMatrix::global_ordinal_type &col, const Scalar &value)
        {
            m.set(row, col, value);
        }

        inline static void add(TpetraMatrix &m, const TpetraMatrix::global_ordinal_type &row, const TpetraMatrix::global_ordinal_type &col, const Scalar &value)
        {
            m.add(row, col, value);
        }

        template<class Tensor>
        inline static void set(Tensor &t, const Scalar value)
        {
            t.set(value);
        }

        template<typename Integer>
        void add_matrix(
            TpetraMatrix &m,
            const std::vector<Integer> &rows,
            const std::vector<Integer> &cols,
            const std::vector<Scalar> &values)
        {
            m.add_matrix(rows, cols, values);
        }

        template<typename Integer>
        void set_matrix(
            TpetraMatrix &m,
            const std::vector<Integer> &rows,
            const std::vector<Integer> &cols,
            const std::vector<Scalar> &values)
        {
            m.set_matrix(rows, cols, values);
        }

        //[host/device locks]
        template<class Tensor>
        static void read_lock(const Tensor &t) {
           const_cast<Tensor &>(t).read_lock();
        }

        template<class Tensor>
        static void read_unlock(const Tensor &t) {
            const_cast<Tensor &>(t).read_unlock();
        }

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

        static void read_and_write_lock(TpetraVector &v)
        {
            v.read_and_write_lock();
        }

        static void read_and_write_unlock(TpetraVector &v)
        {
            v.read_and_write_unlock();
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


        // reductions
        // static Scalar norm2(const TpetraMatrix &m);
        inline static Scalar norm2(const TpetraMatrix &m)
        {
            return m.norm2();
        }

        inline static Scalar norm1(const TpetraMatrix &v)
        {
            assert(false && "IMPLEMENT ME");
            return 0.;
            // return v.norm1();
        }

        static bool is_nan_or_inf(const TpetraMatrix &m)
        {
            assert(false && "IMPLEMENT ME");
            return true;
        }

        inline static Scalar norm_infty(const TpetraMatrix &v)
        {
            assert(false && "IMPLEMENT ME");
            return 0.;
            // return v.norm_infty();
        }

        Scalar reduce(const TpetraMatrix &mat, const Plus &) {
            return mat.sum();
        }

        Scalar reduce(const TpetraVector &vec, const Plus &) {
            return vec.sum();
        }


        Scalar reduce(const TpetraVector &vec, const Max &op) {
            return vec.max();
        }


        Scalar reduce(const TpetraVector &vec, const Min &op) {
            return vec.min();
        }

        void apply_tensor_reduce(TpetraVector &result, const TpetraMatrix &mat, const Plus &, const int dim)
        {
        	if(dim == 1) {
        		TpetraVector vec;
        		vec.values(mat.communicator(), mat.local_size().get(1), mat.size().get(1), 1.);
        		mat.mult(vec, result);
        	} else {
        		TpetraVector vec;
        		vec.values(mat.communicator(), mat.local_size().get(0), mat.size().get(0), 1.);
        		mat.mult_t(vec, result);
        	}
        }

        //blas 1

        inline static void axpy(TpetraVector &y, const Scalar alpha, const TpetraVector &x)
        {
            y.axpy(alpha, x);
        }

        inline static void axpy(TpetraMatrix &y, const Scalar alpha, const TpetraMatrix &x)
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

        inline static void apply_unary(TpetraVector &result, const Minus &, const TpetraVector &v)
        {
            result = v;
            result.scale(-1.);
        }

        template<class Op>
        inline static void apply_unary(TpetraVector &result, const Op &op, const TpetraVector &v)
        {
            result = v;
            result.apply(op);
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
            // assert(!transpose_left);
            //TODO implement transpoe left

            if(transpose_left) {
                left.mult_t(right, result);
            } else {
                left.mult(right, result);
            }
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

        static void apply_binary(TpetraVector &result, const TpetraVector &left, const EMultiplies &, const TpetraVector &right)
        {
            left.e_mul(right, result);
        }

        template<class Op>
        static void apply_binary(TpetraVector &result, const TpetraVector &left, const Op &op, const TpetraVector &right)
        {
            left.apply_binary(op, right, result);
        }

        static void apply_binary(TpetraVector &result, const Reciprocal<Scalar> &reciprocal, const TpetraVector &vec)
        {
            vec.reciprocal(result);

            if(reciprocal.numerator() != 1.) {
                result.scale(reciprocal.numerator());
            }
        }

        // Ac = R*A*P,
        // static void triple_product(
        //     TpetraMatrix &Ac,
        //     const TpetraMatrix &R,
        //     const TpetraMatrix &A,
        //     const TpetraMatrix &P)
        // {
        //     Tpetra::TripleMatrixMultiply::MultiplyRAP(
        //         R,
        //         false, //transposeR
        //         A,
        //         false, //transposeA
        //         P,
        //         false, //transposeP
        //         Ac,
        //         true  //call_FillComplete_on_result
        //     );
        // }

        void diag_scale_left(TpetraMatrix &result, const TpetraVector &diag, const TpetraMatrix &m)
        {
        	result = m;
        	result.implementation().leftScale(diag.implementation());
        }

        static void multiply(
            TpetraMatrix &result,
            bool transpose_left,
            const TpetraMatrix &left,
            bool transpose_right,
            const TpetraMatrix &right)
        {
            left.mult(transpose_left, right, transpose_right, result);
        }

        static void assign_transposed(TpetraMatrix &left, const TpetraMatrix &right)
        {
            right.transpose(left);
        }

        static void diag(TpetraVector &out, const TpetraMatrix &in)
        {
            in.get_diag(out);
        }

        static void diag(TpetraMatrix &out, const TpetraVector &in)
        {
           out.init_diag(in);
        }

        static void diag(TpetraMatrix &out, const TpetraMatrix &in)
        {
            TpetraVector d;
            diag(d, in);
            diag(out, d);
        }

        static bool is_nan_or_inf(const TpetraVector &v)
        {
            return v.is_nan_or_inf();
        }

        template<typename Integer>
        static void set_zero_rows(TpetraMatrix &mat, const std::vector<Integer> &index, const Scalar diag)
        {
            typedef typename TpetraMatrix::global_ordinal_type global_ordinal_type;
            typedef typename TpetraMatrix::Scalar Scalar;

            auto &impl = mat.implementation();
            static const Scalar zero = 0.;

            global_ordinal_type offset = 0;
            Teuchos::ArrayView<const global_ordinal_type> cols;
            Teuchos::ArrayView<const Scalar> values;

            for(auto row : index) {
                if(impl.isGloballyIndexed()) {
                    impl.getGlobalRowView(row, cols, values);
                } else {
                    assert(impl.isLocallyIndexed());
                    auto rr = mat.row_range();
                    impl.getLocalRowView(row - rr.begin(), cols, values);
                    offset = impl.getColMap()->getMinGlobalIndex();
                }

                for(auto c : cols) {
                    const global_ordinal_type col = c + offset;

                    if(col == row) {
                        impl.replaceGlobalValues(row, 1, &diag, &col);
                    } else {
                        impl.replaceGlobalValues(row, 1, &zero, &col);
                    }
                }
            }
        }

        template<typename Integer>
        void select(
            TpetraVector &left,
            const TpetraVector &right,
            const std::vector<Integer> &index)
        {
            right.select(index, left);
        }


        static void read_and_write_lock(TpetraMatrix &t) {
            //IMPLEMENTME
            write_lock(t);
        }

        static void read_and_write_unlock(TpetraMatrix &t){
            //IMPLEMENTME
            write_unlock(t);
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

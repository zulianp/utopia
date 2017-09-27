#ifndef utopia_utopia_BLASBACKEND_CPP
#define utopia_utopia_BLASBACKEND_CPP


#include <vector>
#include <memory>
#include <iostream>
#include <numeric>
#include <algorithm>

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Core.hpp"

#include "utopia_BLASMatrix.hpp"
#include "utopia_BLASTraits.hpp"
#include "utopia_BackendInfo.hpp"
#include "utopia_ScalarBackend.hpp"

#include "utopia_BLASBackend.hpp"


namespace utopia {


    extern "C" {
    // Fortran interface (taken from www.netlib.org/clapack)...
    void dgemm_(const char *transa, const char *transb, const int *m,
                const int *n, const int *k,
                const double *alpha,
                const double *a, const int *lda, const double *b, const int *ldb,
                const double *beta,
                double *c, const int *ldc);

    void dgemv_(const char *trans, const int *m, const int *n,
                const double *alpha, const double *a, const int *lda, const double *x,
                const int *incx, const double *beta, double *c, const int *incy);


    void daxpy_(const int *n,
                const double *alpha,
                const double *x,
                const int *incx,
                double *y,
                const int *incy);

    double dnrm2_(const int *n, const double *x, const int *incx);

    double ddot_(const int *n, const double *sx, const int *incx, const double *sy, const int *incy);

    void dscal_(const int *n,
                const double *sa,
                const double *sx,
                const int *incx);


    #ifdef WITH_MKL

    ///sparse-blas
    void mkl_dcsrgemv_(const char * transa,
                   const int * m,
                   const double * a,
                   const int * ia,
                   const int * ja,
                   const double * x,
                   double * y);

    #endif //WITH_MKL

    }



    void BLASBackend::dscal_wrapper(const int *n, const double *sa, const double *sx, const int *incx)
    {
        dscal_(n, sa, sx, incx); 
    }

    void BLASBackend::daxpy_wrapper(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy)
    {
        daxpy_(n, alpha, x, incx, y, incy); 
    }



    void BLASBackend::assign(Vector &left, Vector &&right)
     {
        using std::move;
        left = move(right);
    }

    void BLASBackend::assign(Vector &left, const Vector &right)
    {
        left = right;
    }


    void BLASBackend::resize(const Size &size, Matrix &mat)
    {
        mat.resize(size.get(0), size.get(1));
    }

    void BLASBackend::resize(const Size &size, Vector &vec)
    {
        vec.resize(size.get(0));
    }

    void BLASBackend::assign(Matrix &left, Matrix &&right)
    {
        using std::move;
        left = move(right);
    }

    void BLASBackend::assign(Matrix &left, const CRSMatrix<Scalar> &right)
    {
        left.resize(right.getRows(), right.getCols());

        std::fill(left.getEntries().begin(), left.getEntries().end(), 0);

        for (SizeType r = 0; r < right.getRows(); ++r) {
            const SizeType begin = right.rowptr()[r];
            const SizeType end   = right.rowptr()[r + 1];

            for (SizeType ind = begin; ind != end; ++ind) {
                const SizeType c = right.colindex()[ind];
                left.set(r, c, right.getEntries()[ind]);
            }
        }
    }

    void BLASBackend::assign(Matrix &left, const Matrix &right)
    {
        left = right;
    }

    void BLASBackend::assign(CCSMatrix<Scalar> &left, const CCSMatrix<Scalar> &right)
     {
        left = right;
    }

    void BLASBackend::assign(CCSMatrix<Scalar> &left, CCSMatrix<Scalar> &&right)
    {
        left = std::move(right);
    }

    void BLASBackend::assign(CRSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right)
    {
        left = right;
    }

    void BLASBackend::assign(CRSMatrix<Scalar> &left, CRSMatrix<Scalar> &&right)
    {
        left = std::move(right);
    }

    void BLASBackend::assign(CCSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right)
     {
        left.resize(right.getRows(), right.getCols());
        CCSMatrix<Scalar>::MapMatrix buffer;
        right.put(buffer, true);
        left.from(buffer);
    }


    void BLASBackend::assignTransposed(Matrix &left, const Matrix &right)
     {
        left.resize(right.getCols(), right.getRows());
        for (SizeType i = 0; i < right.getRows(); ++i) {

            for (SizeType j = 0; j < right.getCols(); ++j) {
                left.set(j, i, right.get(i, j));
            }
        }
    }

    void BLASBackend::assignFromRange(Matrix &left, const Matrix &right, const Range &rowRange, const Range &colRange)
     {

        left.resize(rowRange.extent(), colRange.extent());
        for (SizeType i = 0; i < rowRange.extent(); ++i) {
            for (SizeType j = 0; j < colRange.extent(); ++j) {
                left.set(i, j, right.get(rowRange.begin() + i, colRange.begin() + j));
            }
        }
    }

    void BLASBackend::assignFromRange(Vector &left, const Vector &right, const Range &rowRange, const Range & /*colRange*/)
     {

        left.resize(rowRange.extent());
        for (SizeType i = 0; i < rowRange.extent(); ++i) {
            left[i] = right[rowRange.begin() + i];
        }

    }

    void BLASBackend::select(
                Vector &left,
                const Vector &right,
                const std::vector<SizeType> &index)
    {
        left.resize(index.size());
        for(std::size_t i = 0; i <index.size(); ++i) {
            left[i] = right[index[i]];
        }
    }

    void BLASBackend::select(Matrix &left,
                const Matrix &right,
                const std::vector<SizeType> &row_index,
                const std::vector<SizeType> &col_index)
    {
        

        if(col_index.empty()) {
            left.resize(row_index.size(), right.getCols());
            for(std::size_t i = 0; i < row_index.size(); ++i) {
                for(std::size_t j = 0; j < right.getCols(); ++j) {
                    left.set(i, j, right.get(row_index[i], j));
                }
            }

        } else {
            left.resize(row_index.size(), col_index.size());
            for(std::size_t i = 0; i < row_index.size(); ++i) {
                for(std::size_t j = 0; j < col_index.size(); ++j) {
                    left.set(i, j, right.get(row_index[i], col_index[j]));
                }
            }
        }
    }

    void BLASBackend::assignToRange(Matrix &left, const Matrix &right, const Range &rowRange, const Range &colRange)
     {
        for (unsigned i = 0; i < rowRange.extent(); ++i) {
            for (unsigned j = 0; j < colRange.extent(); ++j) {
                left.set(rowRange.begin() + i, colRange.begin() + j, right.get(i, j));
            }
        }
    }


    void BLASBackend::assignToRange(Matrix &left, const Identity &, const Range &rowRange, const Range &colRange)
     {
        for (unsigned i = 0; i < rowRange.extent(); ++i)
            for (unsigned j = 0; j < colRange.extent(); ++j)
                left.set(rowRange.begin() + i, colRange.begin() + j, i == j);
    }


    void BLASBackend::assignToRange(Vector &left, const Vector &right, const Range &rowRange, const Range & /*colRange*/)
     {
        for (unsigned i = 0; i < rowRange.extent(); ++i)
            left[rowRange.begin() + i] = right[i];
    }

    // inline void BLASBackend::set(Vector &vec, const SizeType index, const BLASBackend::Scalar value)
    //  {
    //     assert(index < vec.size());
    //     vec[index] = value;
    // }

    // inline void BLASBackend::add(Vector &vec, const SizeType index, const BLASBackend::Scalar value)
    //  {
    //     assert(index < vec.size());
    //     vec[index] += value;
    // }


    // inline BLASBackend::Scalar BLASBackend::get(const Vector &vec, const SizeType index)
    //  {
    //     return vec[index];
    // }


    bool BLASBackend::apply(const Vector &left, const Vector &right, const Minus &, Vector &result)
     {
        return zaxpy(-1, right, left, result);
    }

    bool BLASBackend::apply(const Vector &left, const Vector &right, const Plus &, Vector &result)
     {
        return zaxpy(1.0, left, right, result);
    }

    bool BLASBackend::apply(const Matrix &left, const Matrix &right, const Plus &, Matrix &result)
     {
        return zaxpy(1.0, left, right, result);
        //return Matrix();
    }

    bool BLASBackend::apply(const Matrix &left, const Matrix &right, const Multiplies &, Matrix &result)
    {
        return gemm(1, left, right, 0, result);
    }

    bool BLASBackend::apply(const Matrix &left, const Vector &right, const Multiplies &, Vector &result)
     {
        /* Set beta of gemv to 0 so we get y = A * x */
        return gemv(1.0, left, right, 0, result);
    }

    bool BLASBackend::apply(const Matrix &left, const Matrix &right, const Minus &, Matrix &result)
     {
        return zaxpy(-1, right, left, result);
    }

    bool BLASBackend::apply(const Matrix &mat, const Scalar &value, const Multiplies &, Matrix &result)
     {
        return scal(value, mat, result);
    }

    bool BLASBackend::apply(const Matrix &mat, const Scalar &value, const Divides &, Matrix &result)
     {
        return scal(1./value, mat, result);
    }

    bool BLASBackend::apply(const Scalar &value, const CRSMatrix<Scalar> &mat, const Multiplies &, CRSMatrix<Scalar> &result)
     {
        if (&mat != &result) result = mat;
        return scal(value, mat.getEntries(), result.getEntries());
    }


    bool BLASBackend::apply(const Scalar &value, const Matrix &mat, const Multiplies &, Matrix &result)
     {
        return scal(value, mat, result);
    }

    BLASBackend::Scalar BLASBackend::reduce(const Vector &vec, const Plus &)
     {
        using std::accumulate;
        if (vec.empty()) return 0;

        Scalar result = accumulate(vec.begin(), vec.end(), Scalar(0));

        return result;
    }

    BLASBackend::Scalar BLASBackend::reduce(const Matrix &m, const Plus &)
     {
        using std::accumulate;
        if (m.getEntries().empty()) return 0;

//             Compute sum over all elements of the matrix 
        Scalar result = accumulate(m.getEntries().begin(), m.getEntries().end(), Scalar(0));


        return result;
    }



    bool BLASBackend::gemm(const double alpha, const Matrix &left, const Matrix &right, const double beta, Matrix &result)
     {

        const int k = left.getCols();
        assert(k == right.getRows());

        const int m = left.getRows();
        const int n = right.getCols();

        const char Atranspose = 'N';
        const char Btranspose = 'N';

        if (approxeq(beta, 0.0)) {
            result.resize(m, n);
            std::fill(result.getEntries().begin(), result.getEntries().end(), 0);
        }

        dgemm_(&Atranspose, &Btranspose, &m,
               &n, &k,
               &alpha, left.ptr(), &m, right.ptr(), &k,
               &beta, result.ptr(), &m);

        return true;
    }

    bool BLASBackend::gemm(const double alpha, const Matrix &left, const Matrix &right, bool transpose_left,
              bool transpose_right, const double beta, Matrix &result)
               {
        //std::cout << "gemm<double>" << std::endl;

        const int k = transpose_left ? left.getRows() : left.getCols();
        assert(k == (transpose_right ? right.getCols() : right.getRows()));

        const int m = transpose_left ? left.getCols() : left.getRows();
        const int n = transpose_right ? right.getRows() : right.getCols();

        const char Atranspose = transpose_left ? 'T' : 'N';
        const char Btranspose = transpose_right ? 'T' : 'N';

        if (approxeq(beta, 0.0)) {
            result.resize(m, n);
            std::fill(result.getEntries().begin(), result.getEntries().end(), 0);
        }

        dgemm_(&Atranspose, &Btranspose, &m, &n, &k, &alpha, left.ptr(), transpose_left ? &k : &m,
               right.ptr(),
               transpose_right ? &n : &k, &beta, result.ptr(), &m);

        return true;
    }

    bool BLASBackend::gemm(const double scaleFactor, const Vector &left, const Matrix &right, bool transpose_left,
              bool transpose_right, const double beta, Matrix &result)
               {
        assert(transpose_left);
        //std::cout << "gemm<double>" << std::endl;

        const int k = left.size();
        assert(k == (transpose_right ? right.getCols() : right.getRows()));

        const int m = 1;
        const int n = transpose_right ? right.getRows() : right.getCols();

        const char Atranspose = 'T';
        const char Btranspose = transpose_right ? 'T' : 'N';

        if (approxeq(beta, 0.0)) {
            result.resize(m, n);
            std::fill(result.getEntries().begin(), result.getEntries().end(), 0);
        }

        dgemm_(&Atranspose, &Btranspose, &m, &n, &k, &scaleFactor, ptr(left), &k,
               right.ptr(),
               transpose_right ? &n : &k, &beta, result.ptr(), &m);

        return true;
    }

//UNTESTED
      bool BLASBackend::gemm(const double scaleFactor, const Matrix &left, const Vector &right, bool transpose_left,
              bool transpose_right, const double beta, Vector &result)
      {
            assert(!transpose_right);

            const char Atranspose = transpose_left ? 'T' : 'N';

            const int m = left.getRows();
            const int n = left.getCols();

            const int inc = 1;
            const int lda = m;

            result.resize( transpose_left ? left.getCols() : left.getRows() );
            std::fill(result.begin(), result.end(), 0);
      
            //Remember FORTAN counts from 1
            dgemv_(&Atranspose,             //0
                   &m,                      //1
                   &n,                      //2
                   &scaleFactor,            //3
                   ptr(left),               //4
                   &lda,                    //5
                   ptr(right),              //6
                   &inc,                    //7
                   &beta,                   //8
                   ptr(result),             //9
                   &inc);                   //10

            return true;
      }



    bool BLASBackend::apply(const Scalar value, const Vector &vec, const Multiplies &, Vector &result)
     {
        return this->scal(value, vec, result);
    }

    bool BLASBackend::apply(const Vector &vec, const Scalar value, const Multiplies &, Vector &result)
     {
        return this->scal(value, vec, result);
    }


    bool BLASBackend::apply(const Vector &vec, const Scalar value, const Divides &, Vector &result)
     {
        return this->scal(1./value, vec, result);
    }



    bool BLASBackend::gemv(const double alpha, const Matrix &left, const Vector &right, const double beta, Vector &result)
     {

        const char transpose = 'N';

        const int m = left.getRows();
        const int n = left.getCols();

        assert(n == right.size());

        const int incx = 1;
        const int incy = 1;

        /* Check if beta is zero */
        if (approxeq(beta, 0.0)) {
            /* Resize result */
            result.resize(m);
            std::fill(result.begin(), result.end(), 0);
        }

        dgemv_(&transpose, &m, &n,
               &alpha, left.ptr(), &m, &right[0],
               &incx, &beta, &result[0], &incy);

        return true;
    }

    BLASBackend::Scalar BLASBackend::dot(const Vector &left, const Vector &right)
     {
        const int n = left.size();
        assert(n == right.size());

        const int incx = 1;
        const int incy = 1;

        return ddot_(&n, &left[0], &incx, &right[0], &incy);
    }

    BLASBackend::Scalar BLASBackend::norm2(const Vector vector)
     {


        const int n = vector.size();
        const int incx = 1;

        double result = dnrm2_(&n, &vector[0], &incx);

        return result;
    }


     BLASBackend::Scalar BLASBackend::norm_infty(const Vector vector)
      {
        using std::max_element;

        if(vector.empty()) return 0;


       return *max_element(vector.begin(), vector.end());
    }



    bool BLASBackend::size(const Vector &v, Size &size)
     {
        size.setDims(1);
        size.set(0, v.size());
        return true;
    }

    bool BLASBackend::size(const Matrix &m, Size &size)
     {
        size.setDims(2);
        size.set(0, m.getRows());
        size.set(1, m.getCols());
        return true;
    }

    bool BLASBackend::size(const CRSMatrix<Scalar> &m, Size &size)
     {
        size.setDims(2);
        size.set(0, m.getRows());
        size.set(1, m.getCols());
        return true;
    }

    bool BLASBackend::size(const CCSMatrix<Scalar> &m, Size &size)
     {
        size.setDims(2);
        size.set(0, m.getRows());
        size.set(1, m.getCols());
        return true;
    }

    void BLASBackend::build(Matrix &m, const Size &size, const Identity &)
     {
        m.identity(size.get(0), size.get(1));
    }

    void BLASBackend::build(CRSMatrix<Scalar> &m, const Size &size, const Identity &)
     {
        using std::min;

        m.resize(size.get(0), size.get(1));
        const SizeType n = min(m.getRows(), m.getCols());

        m.assemblyBegin();

        for (SizeType i = 0; i != n; ++i) {
            m.set(i, i, 1);
        }

        m.assemblyEnd();
    }

    void BLASBackend::build(CRSMatrix<Scalar> &m, const Size &size, const NNZ<int> &nnz)
    {
        m.initialize(size.get(0), size.get(1), nnz.nnz());
    }

    void BLASBackend::build(CCSMatrix<Scalar> &m, const Size &size, const NNZ<int> &nnz)
    {
        m.initialize(size.get(0), size.get(1), nnz.nnz());
    }

    void BLASBackend::build(CRSMatrix<Scalar> &m, const Size &size, const Zeros & /*values*/)
     {
        m.resize(size.get(0), size.get(1));
    }

    void BLASBackend::build(Matrix &m, const Size &size, const Values<double> &values)
     {
        m.values(size.get(0), size.get(1), values.value());
    }

    void BLASBackend::build(Matrix &m, const Size &size, const Zeros & /*values*/)
     {
        build(m, size, Values<double>(0));
    }

    void BLASBackend::build(Vector &m, const Size &size, const Zeros & /*values*/)
     {
        build(m, size, Values<double>(0));
    }

    void BLASBackend::build(Vector &m, const Size &size, const LocalZeros & /*values*/)
     {
        build(m, size, Values<double>(0));
    }

    void BLASBackend::build(Vector &v, const Size &size, const Values<double> &values)
     {
        v = Vector(size.get(0), values.value());
    }


    bool BLASBackend::mat_diag_shift(Matrix &left, const BLASBackend::Scalar diag_factor)
    {
        using std::min;
        const SizeType n = min(left.getRows(), left.getCols());
        for(SizeType i = 0; i != n; ++i) {
            left.set(i, i, left.get(i, i) + diag_factor);
        }

        return true;
    }


    void BLASBackend::diag(Matrix &left, const Vector &right)
     {
        const SizeType n = right.size();
        left.resize(n, n);
        std::fill(left.getEntries().begin(), left.getEntries().end(), 0);

        for(SizeType i = 0; i != n; ++i) {
            left.set(i, i,right[i]);
        }
    }

    void BLASBackend::diag(Vector &left, const Matrix &right)
    {
        const SizeType n = std::min(right.getRows(), right.getCols());
        left.resize(n);

        for(SizeType i = 0; i != n; ++i) {
            left[i] = right.get(i, i);
        }
    }

    void BLASBackend::diag(Matrix &left, const Matrix &right)
    {
        const SizeType n = std::min(right.getRows(), right.getCols());
        left.resize(n, n);

        std::fill(left.getEntries().begin(), left.getEntries().end(), 0);

        for(SizeType i = 0; i != n; ++i) {
            left.set(i, i, right.get(i, i));
        }
    }

//-------------------------------------------------------------------------------------- //

    bool BLASBackend::apply(const Vector &vec, const Reciprocal<Scalar> &reciprocal, Vector &result)
     {

        result.resize(vec.size());
         for(SizeType i = 0; i < vec.size(); ++i) {
            result[i] = reciprocal.numerator() / vec[i];
          }
          return true;
    }


    bool BLASBackend::outer(const Vector &left, const Vector &right, Matrix &result)
     {
        result.resize(left.size(), right.size());
        for (SizeType i = 0; i < left.size(); ++i) {
            for (SizeType j = 0; j < right.size(); ++j) {
                result.set(i, j, left[i] * right[j]);
            }
        }

        return true;
    }


    bool BLASBackend::apply(const CRSMatrix<double> &left, const Vector &right, const Multiplies &, Vector &result)
     {
        return gemv(1, left, right, result);
    }

    //FIXME use blas 
    bool BLASBackend::transpose(const Matrix &in, Matrix &out)
    {
        out.resize(in.getCols(), in.getRows());
       
        for(SizeType i = 0; i < in.getRows(); ++i) {
            for(SizeType j = 0; j < in.getCols(); ++j) {
                out.set(j, i, in.get(i, j));
            }
        }

        return true;
    }

    inline BLASBackend::Scalar BLASBackend::trace(const Matrix &in) const
    {           

        const SizeType n = std::min(in.getRows(), in.getCols());

        BLASBackend::Scalar ret = 0;
        for(SizeType i = 0; i < n; ++i) {
            ret += in.get(i, i);
        }

        return ret;
    }


#ifdef WITH_MKL

    bool BLASBackend::gemv(const double scaleFactor, const CRSMatrix<double> &left, const Vector &right, Vector &result, const bool transposeMat = false)
    {

        const char transpose = (transposeMat? 'T' : 'N');

        const int m = left.getRows();
        const int n = left.getCols();

        assert(n == right.size());


        result.resize(m);
        std::fill(result.begin(), result.end(), 0);



        ///sparse-blas
        mkl_dcsrgemv_(&transpose,
                       &m,
                       ptr(left),
                       ptr(left.rowptr()),
                       ptr(left.colindex()),
                       ptr(right),
                       ptr(result));

        if(scaleFactor == 1) return true;

       return scal(scaleFactor, result, result);
    }

#else //WITH_MKL

    bool BLASBackend::gemv(const double scaleFactor, const CRSMatrix<double> &left, const Vector &right, Vector &result)
     {
        result.resize(left.getRows());
        for (SizeType r = 0; r != left.getRows(); ++r) {
            result[r] = 0;
            for (SizeType k = left.rowptr()[r]; k != left.rowptr()[r + 1]; ++k) {
                const SizeType c = left.colindex()[k];
                const Scalar val = left.at(k);
                result[r] += val * right[c];
            }

            result[r] *= scaleFactor;
        }

        return true;
    }


#endif //WITH_MKL




}
#endif //utopia_utopia_BLASBACKEND_CPP


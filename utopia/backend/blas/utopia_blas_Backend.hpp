#ifndef utopia_utopia_BLASBACKEND_HPP
#define utopia_utopia_BLASBACKEND_HPP

#include "utopia_blas_Matrix.hpp"
#include "utopia_blas_Traits.hpp"

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Core.hpp"
#include "utopia_BackendInfo.hpp"
#include "utopia_ScalarBackend.hpp"

#include <vector>
#include <memory>
#include <iostream>
#include <numeric>
#include <algorithm>

namespace utopia 
{

    class BLASBackend : public ScalarBackend<double>
    {
        public:
            typedef double Scalar;
            typedef utopia::Matrix<double> Matrix;
            typedef std::vector<double> Vector;
            typedef std::vector<double>::size_type SizeType;

            using ScalarBackend<double>::apply;
            using ScalarBackend<double>::zaxpy;


        void daxpy_wrapper(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
        void dscal_wrapper(const int *n, const double *sa, const double *sx, const int *incx); 


        template<class Tensor>
        inline void readLock(const Tensor &) {}

        template<class Tensor>
        inline void readUnlock(const Tensor &) {}


        template<class Tensor>
        inline void writeLock(const Tensor &) {}

        template<class Tensor>
        inline void writeUnlock(const Tensor &) {}

        inline void writeLock(CRSMatrix<Scalar> &mat) {
            mat.assemblyBegin();
        }

        inline void writeUnlock(CRSMatrix<Scalar> &mat) {
            mat.assemblyEnd();
        }

        inline void writeLock(CCSMatrix<Scalar> &mat) {
            mat.assemblyBegin();
        }

        inline void writeUnlock(CCSMatrix<Scalar> &mat) {
            mat.assemblyEnd();
        }


        template<class Tensor>
        inline void readAndWriteLock(const Tensor &) {}

        template<class Tensor>
        inline void readAndWriteUnlock(const Tensor &) {}

        void assign(Vector &left, Vector &&right);
        void assign(Vector &left, const Vector &right); 

        // template<class V>
        // bool write(const std::string &path, const V &v)
        // {
        //     assert(false);
        //     return false;
        // }
    

        void resize(const Size &size, Matrix &mat); 
        void resize(const Size &size, Vector &vec); 
        void assign(Matrix &left, Matrix &&right) ; 
        void assign(Matrix &left, const CRSMatrix<Scalar> &right);  
        void assign(Matrix &left, const Matrix &right);
        void assign(CCSMatrix<Scalar> &left, const CCSMatrix<Scalar> &right); 
        void assign(CCSMatrix<Scalar> &left, CCSMatrix<Scalar> &&right);
        void assign(CRSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right);
        void assign(CRSMatrix<Scalar> &left, CRSMatrix<Scalar> &&right);
        void assign(CCSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right);
    

        template<class Tensor>
        void build(Tensor &t, const Size &s, const Resize &)
        {
            resize(s, t);
        }


        inline static Range range(const Vector &v)
        {
            return Range(0, v.size());
        }

        

        template<class MatrixT>
        inline static Range rowRange(const MatrixT &m) {
            return Range(0, m.getRows());
        }

        template<class MatrixT>
        inline static Range colRange(const MatrixT &m) {
            return Range(0, m.getCols());
        }

        void assignTransposed(Matrix &left, const Matrix &right);
        void assignFromRange(Matrix &left, const Matrix &right, const Range &rowRange, const Range &colRange);
        void assignFromRange(Vector &left, const Vector &right, const Range &rowRange, const Range & /*colRange*/);
        void select(Vector &left,
                    const Vector &right,
                    const std::vector<SizeType> &index);

        void select(Matrix &left,
                    const Matrix &right,
                    const std::vector<SizeType> &row_index,
                    const std::vector<SizeType> &col_index);

        void assignToRange(Matrix &left, const Matrix &right, const Range &rowRange, const Range &colRange);
        void assignToRange(Matrix &left, const Identity &, const Range &rowRange, const Range &colRange);
        void assignToRange(Vector &left, const Vector &right, const Range &rowRange, const Range & /*colRange*/);
        

        inline void set(Vector &vec, const SizeType index, const Scalar value)
         {
            assert(index < vec.size());
            vec[index] = value;
        }

        inline void add(Vector &vec, const SizeType index, const Scalar value)
         {
            assert(index < vec.size());
            vec[index] += value;
        }


        inline Scalar get(const Vector &vec, const SizeType index)
         {
            return vec[index];
        }

        

        template<class MatrixT>
        inline Scalar get(const MatrixT &mat, const SizeType row, const SizeType col) {
        //            return mat.get(row, col);

            assert(row < mat.getRows() && col < mat.getCols());
            return mat.get(row, col);
        }

        template<class MatrixT>
        inline void set(MatrixT &mat, const SizeType row, const SizeType col, const Scalar value) {
            assert(row < mat.getRows() && col < mat.getCols());
            mat.set(row, col, value);
        }

        template<class MatrixT>
        inline void add(MatrixT &mat, const SizeType row, const SizeType col, const Scalar value) {
            assert(row < mat.getRows() && col < mat.getCols());
            mat.set(row, col, mat.get(row, col) + value);
        }

        bool apply(const Vector &left, const Vector &right, const Minus &, Vector &result);
        bool apply(const Vector &left, const Vector &right, const Plus &, Vector &result);
        bool apply(const Matrix &left, const Matrix &right, const Plus &, Matrix &result);
        bool apply(const Matrix &left, const Matrix &right, const Multiplies &, Matrix &result); 
        bool apply(const Matrix &left, const Vector &right, const Multiplies &, Vector &result);
        bool apply(const Matrix &left, const Matrix &right, const Minus &, Matrix &result);
        bool apply(const Matrix &mat, const Scalar &value, const Multiplies &, Matrix &result);
        bool apply(const Matrix &mat, const Scalar &value, const Divides &, Matrix &result);
        bool apply(const Scalar &value, const CRSMatrix<Scalar> &mat, const Multiplies &, CRSMatrix<Scalar> &result);
        bool apply(const Scalar &value, const Matrix &mat, const Multiplies &, Matrix &result);
        
        Scalar reduce(const Vector &vec, const Plus &);
        Scalar reduce(const Matrix &m, const Plus &);
         
        ///////////////////////////////////////////////////////////  check ////////////////////////////
        template<class Left, class Right, class Result>
        bool zaxpy(const double scaleFactor, const Left &left, Right &&right, Result &result) 
        {
            using std::forward;

            assert(left.size() == right.size());

            result = forward<Right>(right);
            const int isize = result.size();
            const int incx = 1;
            const int incy = 1;
            //call to blas
            // daxpy_(&isize, &scaleFactor, ptr(left), &incx, ptr(result), &incy);
            daxpy_wrapper(&isize, &scaleFactor, ptr(left), &incx, ptr(result), &incy);

            return true;
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////


        bool gemm(const double alpha, const Matrix &left, const Matrix &right, const double beta, Matrix &result);
        bool gemm(const double alpha, const Matrix &left, const Matrix &right, bool transpose_left, bool transpose_right, const double beta, Matrix &result);
        bool gemm(const double scaleFactor, const Vector &left, const Matrix &right, bool transpose_left, bool transpose_right, const double beta, Matrix &result);
        

    //UNTESTED
        bool gemm(const double scaleFactor, const Matrix &left, const Vector &right, bool transpose_left, bool transpose_right, const double beta, Vector &result);

        template<class Left, class Right, class Result>
        bool gem(const double scaleFactor, const Left &left, const Right &right, bool transpose_left,
                  bool transpose_right, const double beta, Result &result)
        {
            return gemm(scaleFactor, left, right, transpose_left, transpose_right, beta, result);
        }


        inline static Scalar *ptr(Matrix &mat)
         {
            return mat.ptr();
        }

        inline static Scalar *ptr(Vector &v)
         {
            return &v[0];
        }

        inline static const Scalar *ptr(const Matrix &mat)
         {
            return mat.ptr();
        }

        inline static const Scalar *ptr(const Vector &v)
         {
            return &v[0];
        }

        inline static const Scalar *ptr(const CRSMatrix<Scalar> &mat)
         {
            return &mat.getEntries()[0];
        }


        template<typename T>
        inline static const T *ptr(const std::vector<T> &v) {
            return &v[0];
        }

        bool apply(const Scalar value, const Vector &vec, const Multiplies &, Vector &result);
        bool apply(const Vector &vec, const Scalar value, const Multiplies &, Vector &result);
        bool apply(const Vector &vec, const Scalar value, const Divides &, Vector &result);
        

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        template<class VectorT>
        bool scal(const double scaleFactor, VectorT &&v, Vector &result) 
        {

            const int n = v.size();
            const int incx = 1;

            result = std::forward<VectorT>(v);
            assert(!result.empty());

            // dscal_(&n, &scaleFactor, &result[0], &incx);
            dscal_wrapper(&n, &scaleFactor, &result[0], &incx);

            return true;

        }

        template<class MatrixT>
        bool scal(const double scaleFactor, MatrixT &&m, Matrix &result) 
        {

            const int n = m.getRows() * m.getCols();
            const int incx = 1;
            result = std::forward<MatrixT>(m);
            dscal_wrapper(&n, &scaleFactor, result.ptr(), &incx);

            return true;
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        

        bool gemv(const double alpha, const Matrix &left, const Vector &right, const double beta, Vector &result);
        Scalar dot(const Vector &left, const Vector &right);
        Scalar norm2(const Vector vector);
        Scalar norm_infty(const Vector vector);
        

        template<class Tensor>
        bool local_size(const Tensor &t, Size &size)
        {
            return this->size(t, size);
        }


        bool size(const Vector &v, Size &size);
        bool size(const Matrix &m, Size &size);
        bool size(const CRSMatrix<Scalar> &m, Size &size);
        bool size(const CCSMatrix<Scalar> &m, Size &size);

        void build(Matrix &m, const Size &size, const Identity &);
        void build(CRSMatrix<Scalar> &m, const Size &size, const Identity &);
        void build(CRSMatrix<Scalar> &m, const Size &size, const NNZ<int> &nnz);
        void build(CCSMatrix<Scalar> &m, const Size &size, const NNZ<int> &nnz);
        void build(CRSMatrix<Scalar> &m, const Size &size, const Zeros & /*values*/);
        void build(Matrix &m, const Size &size, const Values<double> &values);
        void build(Matrix &m, const Size &size, const Zeros & /*values*/);
        void build(Vector &m, const Size &size, const Zeros & /*values*/);
        void build(Vector &m, const Size &size, const LocalZeros & /*values*/);
        void build(Vector &v, const Size &size, const Values<double> &values);



        bool mat_diag_shift(Matrix &left, const Scalar diag_factor);
        void diag(Matrix &left, const Vector &right);
        void diag(Vector &left, const Matrix &right);
        void diag(Matrix &left, const Matrix &right);
    

//------------------------------------------TO BE DONE -------------------------------- //
        template<class Tensor>
         bool monitor(const long &, Tensor &)
        {
            std::cout<<"BLAS_backend:: monitor needs to be implemented  \n";     
            return true; 
        }
        
        template<class Tensor>
        bool set_zero_rows(Tensor & /*Mat_A */, const std::vector<int> & /*index*/)
        {
            std::cout<<"BLAS_backend:: set_zero_rows needs to be implemented  \n";     
            return true; 
        }

        template<class MatTensor, class VecTensor>
        bool apply_BC_to_system(MatTensor & /*A*/, VecTensor& /*x*/, VecTensor& /*rhs*/, const std::vector<int> & /*index*/)
        {
            std::cout<<"BLAS_backend:: apply_BC_to_system needs to be implemented  \n";     
            return true; 
        }


        // TODO:: 
        template<typename Matrix>
        bool triple_product_PtAP(const Matrix & /*M*/, const Matrix & /*I*/, Matrix & /*result*/)
        {
            std::cout<<"BLAS:: triple_product_PtAP not implemented ...  \n";     
            return true; 
        }
//-------------------------------------------------------------------------------------- //

        
        static bool diag_scale_right(const Matrix &m, const Vector &diag, Matrix &result)
        {
            result.resize(m.getRows(), m.getCols());

            const SizeType n = diag.size();
            ASSERT(n == m.getCols() && "sizes are not compatible");

            for(SizeType i = 0; i < m.getRows(); ++i) {
                for(SizeType j = 0; j < m.getCols(); ++j) {
                    result.set(i, j, m.get(i, j) * diag[j]);
                }
            }

            return true;
        }

        static bool diag_scale_left(const Vector &diag, const Matrix &m, Matrix &result)
        {

            result.resize(m.getRows(), m.getCols());

            const SizeType n = diag.size();
            ASSERT(n == m.getRows() && "sizes are not compatible");

            for(SizeType i = 0; i < m.getRows(); ++i) {
                for(SizeType j = 0; j < m.getCols(); ++j) {
                    result.set(i, j, m.get(i, j) * diag[i]);
                }
            }

            return true;
        }



        template<typename Ordinal>
        inline void set(Vector &v, const std::vector<Ordinal> &indices, const std::vector<Scalar> &values) {
            assert(indices.size() == values.size());

            for (typename std::vector<Ordinal>::size_type i = 0; i < indices.size(); ++i)
                set(v, indices[i], values[i]);

        }

        template<typename Ordinal>
        inline void set(Matrix &m, const std::vector<Ordinal> &rows, const std::vector<Ordinal> &columns,
                        const std::vector<Scalar> &values) {
            assert(rows.size() == columns.size() && rows.size() == values.size());

            for (typename std::vector<Ordinal>::size_type i = 0; i < rows.size(); ++i)
                set(m, rows[i], columns[i], values[i]);
        }

        template<typename Ordinal>
        inline void set(CRSMatrix<Scalar> &m, const std::vector<Ordinal> &rows, const std::vector<Ordinal> &columns,
                        const std::vector<Scalar> &values) {

            for (typename std::vector<Ordinal>::size_type i = 0; i < rows.size(); ++i) {
                set(m, rows[i], columns[i], values[i]);
            }
        }

        template<typename T, class Comparator>
        bool compare(const std::vector<T> &left, const std::vector<T> &right, const Comparator &comp) {
            //DEBUG
            assert(left.size() == right.size());

            //RELEASE
            if (left.size() != right.size()) return false;

            bool result = true;
            for (SizeType i = 0; i < left.size(); ++i) {
                result &= comp(left[i], right[i]);
            }

            return result;
        }

        template<class Comparator>
        bool compare(const Matrix &left, const Matrix &right, const Comparator &comp) {
            //DEBUG
            assert(left.getRows() == right.getRows());
            assert(left.getCols() == right.getCols());

            //RELEASE
            if (left.getRows() != right.getRows() || left.getCols() != right.getCols()) return false;

            return compare(left.getEntries(), right.getEntries(), comp);
        }


        template<class Comparator>
        bool compare(const CRSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right, const Comparator &comp) {
            //DEBUG
            assert(left.getRows() == right.getRows());
            assert(left.getCols() == right.getCols());

            //RELEASE
            if (left.getRows() != right.getRows() || left.getCols() != right.getCols()) return false;


            for (SizeType r = 0; r != left.getRows(); ++r) {
                const SizeType boundaryl = left.rowptr()[r + 1];
                const SizeType boundaryr = right.rowptr()[r + 1];

                for (SizeType kl = left.rowptr()[r], kr = right.rowptr()[r];
                     kl != boundaryl && kr != boundaryr;
                    /* increment inside */
                        )
                {
                    const SizeType cl = left.colindex()[kl];
                    const SizeType cr = right.colindex()[kr];
                    const Scalar vl = left.at(kl);
                    const Scalar vr = right.at(kr);

                    if(cl == cr) {
                        if (!comp(vl, vr)) {
                            return false;
                        } else {
                            ++kl;
                            ++kr;
                        }
                    }
                    else if(cl < cr)
                    {
                        if(!comp(vl, Scalar(0))) {
                            return false;
                        }
                        ++kl;
                    }
                    else if(cl > cr) {
                        if(!comp(vr, Scalar(0))) {
                            return false;
                        }
                        ++kr;
                    }

                    if(kl == boundaryl && kr != boundaryr) {
                        for(; kl != boundaryl; ++kl) {
                            if(!comp(left.at(kl), Scalar(0))) {
                                return false;
                            }
                        }

                    } else if(kl != boundaryl && kr == boundaryr) {
                        for(; kr != boundaryr; ++kr) {
                            if(!comp(right.at(kr), Scalar(0))) {
                                return false;
                            }
                        }
                    }
                }
            }

           return true;
        }

        template<class Operation>
        bool apply(const Vector &vec, const Operation &, Vector &result) {


            result.resize(vec.size());
            for(SizeType i = 0; i < vec.size(); ++i) {
                result[i] = Operation::template apply<Scalar>(vec[i]);
            }

            return true;
        }

        bool apply(const Vector &vec, const Minus &, Vector &result) {
            result.resize(vec.size());
            for(SizeType i = 0; i < vec.size(); ++i) {
                result[i] = -vec[i];
            }

            return true;
        }



        template<class Operation>
        bool apply(const Vector &left, const Vector &right, const Operation &op, Vector &result) {


            result.resize(left.size());
            return apply(left, right, result, op);
        }


        bool apply(const Vector &vec, const Reciprocal<Scalar> &reciprocal, Vector &result);


        template<class Operation>
        bool apply(const Vector &left, const Vector &right, Vector &result, const Operation &op) {
            for (SizeType i = 0; i < left.size(); ++i) {
                result[i] = op.template apply<Scalar>(left[i], right[i]);
            }
            return true;
        }


        bool handle(int err) const {
            ASSERT(err == 0);
            return true;
        }

        bool outer(const Vector &left, const Vector &right, Matrix &result);
        bool apply(const CRSMatrix<double> &left, const Vector &right, const Multiplies &, Vector &result);

        //FIXME use blas 
        bool transpose(const Matrix &in, Matrix &out);
        inline Scalar trace(const Matrix &in) const;



#ifdef WITH_MKL

        bool gemv(const double scaleFactor, const CRSMatrix<double> &left, const Vector &right, Vector &result, const bool transposeMat = false);


#else //WITH_MKL

        bool gemv(const double scaleFactor, const CRSMatrix<double> &left, const Vector &right, Vector &result);


#endif //WITH_MKL

    };


    template<typename T>
    class Backend<T, BLAS> : public BLASBackend {
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
            info_.set_name("blas");
        }
    };
}
#endif //utopia_utopia_BLASBACKEND_HPP
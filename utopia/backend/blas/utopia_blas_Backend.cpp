// #ifndef utopia_utopia_BLASBACKEND_CPP
// #define utopia_utopia_BLASBACKEND_CPP

// #include "utopia_blas_Matrix.hpp"
// #include "utopia_blas_Traits.hpp"
// #include "utopia_blas_Backend.hpp"

// #include "utopia_ForwardDeclarations.hpp"
// #include "utopia_Wrapper.hpp"
// #include "utopia_Core.hpp"
// #include "utopia_BackendInfo.hpp"
// #include "utopia_ScalarBackend.hpp"

// #include <vector>
// #include <memory>
// #include <iostream>
// #include <numeric>
// #include <algorithm>
// #include <fstream>
// #include <string>

// namespace utopia {


//     extern "C" {
//         // Fortran interface (taken from www.netlib.org/clapack)...
//         void dgemm_(const char *transa, const char *transb, const int *m,
//                     const int *n, const int *k,
//                     const double *alpha,
//                     const double *a, const int *lda, const double *b, const int *ldb,
//                     const double *beta,
//                     double *c, const int *ldc);

//         void dgemv_(const char *trans, const int *m, const int *n,
//                     const double *alpha, const double *a, const int *lda, const double *x,
//                     const int *incx, const double *beta, double *c, const int *incy);


//         void daxpy_(const int *n,
//                     const double *alpha,
//                     const double *x,
//                     const int *incx,
//                     double *y,
//                     const int *incy);

//         double dnrm2_(const int *n, const double *x, const int *incx);
//         double dasum_(const int *n, const double *x, const int *incx);

//         double ddot_(const int *n, const double *sx, const int *incx, const double *sy, const int *incy);

//         void dscal_(const int *n,
//                     const double *sa,
//                     const double *sx,
//                     const int *incx);


// #ifdef WITH_MKL

//         ///sparse-blas
//         void mkl_dcsrgemv_(const char * transa,
//                            const int * m,
//                            const double * a,
//                            const int * ia,
//                            const int * ja,
//                            const double * x,
//                            double * y);

// #endif //WITH_MKL

//     }

//     void BLASBackend::dscal_wrapper(const int *n, const double *sa, const double *sx, const int *incx)
//     {
//         dscal_(n, sa, sx, incx);
//     }

//     void BLASBackend::daxpy_wrapper(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy)
//     {
//         daxpy_(n, alpha, x, incx, y, incy);
//     }


//     void BLASBackend::scale(Vector &result, const Scalar scale_factor)
//     {
//         const int n = result.size();
//         const int incx = 1;
//         assert(!result.empty());
//         dscal_wrapper(&n, &scale_factor, &result[0], &incx);
//     }

//     void BLASBackend::scale(Matrix &result, const Scalar scale_factor)
//     {
//         const int n = result.rows() * result.cols();
//         const int incx = 1;
//         dscal_wrapper(&n, &scale_factor, result.ptr(), &incx);
//     }

//     void BLASBackend::assign(Vector &left, Vector &&right)
//     {
//         using std::move;
//         left = move(right);
//     }

//     void BLASBackend::assign(Vector &left, const Vector &right)
//     {
//         left = right;
//     }

//     void BLASBackend::resize(Matrix &mat, const Size &size)
//     {
//         mat.resize(size.get(0), size.get(1));
//     }

//     void BLASBackend::resize(Vector &vec, const Size &size)
//     {
//         vec.resize(size.get(0));
//     }

//     void BLASBackend::assign(Matrix &left, Matrix &&right)
//     {
//         using std::move;
//         left = move(right);
//     }

//     void BLASBackend::assign(Matrix &left, const CRSMatrix<Scalar> &right)
//     {
//         left.resize(right.rows(), right.cols());

//         std::fill(left.entries().begin(), left.entries().end(), 0);

//         for (decltype(right.rows()) r = 0; r < right.rows(); ++r) {
//             const SizeType begin = right.rowptr()[r];
//             const SizeType end   = right.rowptr()[r + 1];

//             for (SizeType ind = begin; ind != end; ++ind) {
//                 const SizeType c = right.colindex()[ind];
//                 left.set(r, c, right.entries()[ind]);
//             }
//         }
//     }

//     void BLASBackend::assign(Matrix &left, const Matrix &right)
//     {
//         left = right;
//     }

//     void BLASBackend::assign(CCSMatrix<Scalar> &left, const CCSMatrix<Scalar> &right)
//     {
//         left = right;
//     }

//     void BLASBackend::assign(CCSMatrix<Scalar> &left, CCSMatrix<Scalar> &&right)
//     {
//         left = std::move(right);
//     }

//     void BLASBackend::assign(CRSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right)
//     {
//         left = right;
//     }

//     void BLASBackend::assign(CRSMatrix<Scalar> &left, CRSMatrix<Scalar> &&right)
//     {
//         left = std::move(right);
//     }

//     void BLASBackend::assign(CCSMatrix<Scalar> &left, const CRSMatrix<Scalar> &right)
//     {
//         left.resize(right.rows(), right.cols());
//         CCSMatrix<Scalar>::MapMatrix buffer;
//         right.put(buffer, true);
//         left.from(buffer);
//     }

//     void BLASBackend::assign_from_range(Matrix &left, const Matrix &right, const Range &rowRange, const Range &colRange)
//     {

//         left.resize(rowRange.extent(), colRange.extent());

//         const SizeType r_ext = rowRange.extent();
//         const SizeType c_ext = colRange.extent();

//         for (SizeType i = 0; i < r_ext; ++i) {
//             for (SizeType j = 0; j < c_ext; ++j) {
//                 left.set(i, j, right.get(rowRange.begin() + i, colRange.begin() + j));
//             }
//         }
//     }

//     void BLASBackend::assign_from_range(Vector &left, const Vector &right, const Range &rowRange, const Range & /*colRange*/)
//     {
//         const SizeType r_ext = rowRange.extent();

//         left.resize(r_ext);
//         for (SizeType i = 0; i < r_ext; ++i) {
//             left[i] = right[rowRange.begin() + i];
//         }

//     }

//     void BLASBackend::select(Vector &left,
//                              const Vector &right,
//                              const std::vector<SizeType> &index)
//     {
//         left.resize(index.size());
//         for(std::size_t i = 0; i <index.size(); ++i) {
//             left[i] = right[index[i]];
//         }
//     }

//     void BLASBackend::select(Matrix &left,
//                              const Matrix &right,
//                              const std::vector<SizeType> &row_index,
//                              const std::vector<SizeType> &col_index)
//     {


//         if(col_index.empty()) {
//             left.resize(row_index.size(), right.cols());
//             for(std::size_t i = 0; i < row_index.size(); ++i) {
//                 for(std::size_t j = 0; j < right.cols(); ++j) {
//                     left.set(i, j, right.get(row_index[i], j));
//                 }
//             }

//         } else {
//             left.resize(row_index.size(), col_index.size());
//             for(std::size_t i = 0; i < row_index.size(); ++i) {
//                 for(std::size_t j = 0; j < col_index.size(); ++j) {
//                     left.set(i, j, right.get(row_index[i], col_index[j]));
//                 }
//             }
//         }
//     }

//     void BLASBackend::assign_to_range(Matrix &left, const Matrix &right, const Range &rowRange, const Range &colRange)
//     {
//         for (unsigned i = 0; i < rowRange.extent(); ++i) {
//             for (unsigned j = 0; j < colRange.extent(); ++j) {
//                 left.set(rowRange.begin() + i, colRange.begin() + j, right.get(i, j));
//             }
//         }
//     }


//     void BLASBackend::assign_to_range(Matrix &left, const Identity &, const Range &rowRange, const Range &colRange)
//     {
//         for (unsigned i = 0; i < rowRange.extent(); ++i) {
//             for (unsigned j = 0; j < colRange.extent(); ++j) {
//                 left.set(rowRange.begin() + i, colRange.begin() + j, i == j);
//             }
//         }
//     }

//     void BLASBackend::assign_to_range(Vector &left, const Vector &right, const Range &rowRange, const Range & /*colRange*/)
//     {
//         for (unsigned i = 0; i < rowRange.extent(); ++i) {
//             left[rowRange.begin() + i] = right[i];
//         }
//     }

//     BLASBackend::Scalar BLASBackend::reduce(const Vector &vec, const Plus &)
//     {
//         using std::accumulate;
//         if (vec.empty()) return 0.;
//         return accumulate(vec.begin(), vec.end(), Scalar(0));
//     }

//     BLASBackend::Scalar BLASBackend::reduce(const Matrix &m, const Plus &)
//     {
//         using std::accumulate;
//         if (m.empty()) return 0.;
//         return accumulate(m.entries().begin(), m.entries().end(), Scalar(0));
//     }

//     void BLASBackend::gemm(Matrix &C,
//                            const BLASBackend::Scalar beta,
//                            const BLASBackend::Scalar &alpha,
//                            const bool transpose_A,
//                            const Matrix &A,
//                            const bool transpose_B,
//                            const Matrix &B) {

//         const int k = transpose_A ? A.rows() : A.cols();
//         assert(k == (transpose_B ? B.cols() : B.rows()));

//         const int m = transpose_A ? A.cols() : A.rows();
//         const int n = transpose_B ? B.rows() : B.cols();

//         const char t_A_flag = transpose_A ? 'T' : 'N';
//         const char t_B_flag = transpose_B ? 'T' : 'N';

//         if (C.empty() || approxeq(beta, 0.0)) {
//             C.resize(m, n);
//             std::fill(C.entries().begin(), C.entries().end(), 0);
//         }

//         dgemm_(&t_A_flag, &t_B_flag, &m, &n, &k, &alpha, A.ptr(), transpose_A ? &k : &m,
//                B.ptr(),
//                transpose_B ? &n : &k, &beta, C.ptr(), &m);
//     }


//     void BLASBackend::gemm(Matrix &C,
//                            const Scalar beta,
//                            const Scalar &alpha,
//                            const bool transpose_A,
//                            const Vector &A,
//                            const bool transpose_B,
//                            const Matrix &B)
//     {

//         const int k = transpose_A ? A.size() : 1;
//         assert(k == (transpose_B ? B.cols() : B.rows()));

//         const int m = transpose_A ? 1 : A.size();
//         const int n = transpose_B ? B.rows() : B.cols();

//         const char t_A_flag = transpose_A ? 'T' : 'N';
//         const char t_B_flag = transpose_B ? 'T' : 'N';

//         if (C.empty() || approxeq(beta, 0.0)) {
//             C.resize(m, n);
//             std::fill(C.entries().begin(), C.entries().end(), 0);
//         }

//         dgemm_(&t_A_flag, &t_B_flag, &m, &n, &k, &alpha, ptr(A), transpose_A ? &k : &m,
//                B.ptr(),
//                transpose_B ? &n : &k, &beta, C.ptr(), &m);
//     }

//     void BLASBackend::gemv(Vector &y,
//                            const Scalar beta,
//                            const Scalar &alpha,
//                            const bool transpose_A,
//                            const Matrix &A,
//                            const Vector &x)
//     {
//         const char t_A_flag = transpose_A ? 'T' : 'N';

//         const int m = A.rows();
//         const int n = A.cols();

//         const int y_size = transpose_A ? A.cols() : A.rows();

//         const int inc = 1;
//         const int lda = m;

//         if(y.empty() || y_size != int(y.size())) {
//             y.resize(y_size);
//             std::fill(y.begin(), y.end(), 0);
//         } else if(approxeq(alpha, 0.0)) {
//             std::fill(y.begin(), y.end(), 0);
//         }

//         dgemv_(&t_A_flag,               //0
//                &m,                      //1
//                &n,                      //2
//                &alpha,            //3
//                ptr(A),               //4
//                &lda,                    //5
//                ptr(x),              //6
//                &inc,                    //7
//                &beta,                   //8
//                ptr(y),                  //9
//                &inc);                   //10

//     }

//     void BLASBackend::gemv(Vector &y,
//                            const Scalar beta,
//                            const Scalar &alpha,
//                            const bool transpose_A,
//                            const CRSMatrix<Scalar> &A,
//                            const Vector &x)
//     {
//         int y_size = transpose_A? A.rows() : A.cols();

//         if(y.empty() || int(y.size()) != y_size) {
//             y.resize(y_size);
//             std::fill(y.begin(), y.end(), 0);
//         }

//         const SizeType A_rows = A.rows();

//         if(transpose_A) {
//             for (SizeType r = 0; r != A_rows; ++r) {
//                 BLASBackend::Scalar x_r = x[r];

//                 for (auto k = A.rowptr()[r]; k != A.rowptr()[r + 1]; ++k) {
//                     const SizeType c = A.colindex()[k];
//                     const Scalar A_rk = A.at(k);
//                     y[c] += alpha * A_rk * x_r;
//                 }
//             }

//         } else {
//             for (SizeType r = 0; r != A_rows; ++r) {
//                 BLASBackend::Scalar A_x = 0.;
//                 for (auto k = A.rowptr()[r]; k != A.rowptr()[r + 1]; ++k) {
//                     const SizeType c = A.colindex()[k];
//                     const Scalar A_rk = A.at(k);
//                     A_x += A_rk * x[c];
//                 }

//                 y[r] += alpha * A_x;
//             }
//         }
//     }

//     void BLASBackend::axpy(Vector &y, const Scalar &alpha, const Vector &x)
//     {
//         if(y.size() != x.size()) {
//             y.resize(x.size());
//             std::fill(y.begin(), y.end(), 0.);
//         }

//         int n = x.size();
//         int incx = 1;
//         int incy = 1;

//         daxpy_wrapper(&n, &alpha, ptr(x), &incx, ptr(y), &incy);
//     }

//     void BLASBackend::axpy(Matrix &y, const Scalar &alpha, const Matrix &x)
//     {
//         if(y.rows() != x.rows() || y.cols() != x.cols()) {
//             y.resize(x.rows(), x.cols());
//             std::fill(y.entries().begin(), y.entries().end(), 0.);
//         }

//         axpy(y.entries(), alpha, x.entries());
//     }

//     BLASBackend::Scalar BLASBackend::dot(const Vector &left, const Vector &right)
//     {
//         const int n = left.size();
//         assert(n == right.size());

//         const int incx = 1;
//         const int incy = 1;

//         return ddot_(&n, &left[0], &incx, &right[0], &incy);
//     }

//     BLASBackend::Scalar BLASBackend::dot(const Matrix &left, const Matrix &right)
//     {
//         const int n = left.n_elements();
//         assert(n == right.n_elements());

//         const int incx = 1;
//         const int incy = 1;

//         return ddot_(&n, &left.entries()[0], &incx, &right.entries()[0], &incy);
//     }

//     BLASBackend::Scalar BLASBackend::norm2(const Vector &vector)
//     {
//         const int n = vector.size();
//         const int incx = 1;
//         return dnrm2_(&n, &vector[0], &incx);
//     }

//     BLASBackend::Scalar BLASBackend::norm1(const Vector &vector)
//     {
//         const int n = vector.size();
//         const int incx = 1;
//         return dasum_(&n, &vector[0], &incx);
//     }

//     BLASBackend::Scalar BLASBackend::norm2(const Matrix &mat)
//     {
//         const int n = mat.entries().size();
//         const int incx = 1;
//         return dnrm2_(&n, &mat.entries()[0], &incx);
//     }

//     BLASBackend::Scalar BLASBackend::norm2(const CRSMatrix<Scalar> &mat)
//     {
//         const int n = mat.entries().size();
//         const int incx = 1;
//         return dnrm2_(&n, &mat.entries()[0], &incx);
//     }

//     BLASBackend::Scalar BLASBackend::norm2(const CCSMatrix<Scalar> &mat)
//     {
//         const int n = mat.entries().size();
//         const int incx = 1;
//         return dnrm2_(&n, &mat.entries()[0], &incx);
//     }

//     BLASBackend::Scalar BLASBackend::norm_infty(const Vector &vector)
//     {
//         using std::max_element;
//         if(vector.empty()) return 0.;
//         return *max_element(vector.begin(), vector.end());
//     }

//     void BLASBackend::size(const Vector &v, Size &size)
//     {
//         size.set_dims(1);
//         size.set(0, v.size());
//     }

//     void BLASBackend::size(const Matrix &m, Size &size)
//     {
//         size.set_dims(2);
//         size.set(0, m.rows());
//         size.set(1, m.cols());
//     }

//     void BLASBackend::size(const CRSMatrix<Scalar> &m, Size &size)
//     {
//         size.set_dims(2);
//         size.set(0, m.rows());
//         size.set(1, m.cols());
//     }

//     void BLASBackend::size(const CCSMatrix<Scalar> &m, Size &size)
//     {
//         size.set_dims(2);
//         size.set(0, m.rows());
//         size.set(1, m.cols());
//     }

//     void BLASBackend::build(Matrix &m, const Size &size, const Identity &)
//     {
//         m.identity(size.get(0), size.get(1));
//     }

//     void BLASBackend::build(Matrix &m, const Size &size, const LocalIdentity &)
//     {
//         m.identity(size.get(0), size.get(1));
//     }


//     void BLASBackend::build(CRSMatrix<Scalar> &m, const Size &size, const Identity &)
//     {
//         using std::min;

//         m.resize(size.get(0), size.get(1));
//         const SizeType n = min(m.rows(), m.cols());

//         m.assembly_begin();

//         for (SizeType i = 0; i != n; ++i) {
//             m.set(i, i, 1);
//         }

//         m.assembly_end();
//     }

//     void BLASBackend::build(CRSMatrix<Scalar> &m, const Size &size, const LocalIdentity &)
//     {
//         using std::min;

//         m.resize(size.get(0), size.get(1));
//         const SizeType n = min(m.rows(), m.cols());

//         m.assembly_begin();

//         for (SizeType i = 0; i != n; ++i) {
//             m.set(i, i, 1);
//         }

//         m.assembly_end();
//     }

//     void BLASBackend::build(CRSMatrix<Scalar> &m, const Size &size, const NNZ<int> &nnz)
//     {
//         m.initialize(size.get(0), size.get(1), nnz.nnz());
//     }

//     void BLASBackend::build(CCSMatrix<Scalar> &m, const Size &size, const NNZ<int> &nnz)
//     {
//         m.initialize(size.get(0), size.get(1), nnz.nnz());
//     }

//     void BLASBackend::build(CRSMatrix<Scalar> &m, const Size &size, const Zeros & /*values*/)
//     {
//         m.resize(size.get(0), size.get(1));
//     }

//     void BLASBackend::build(Matrix &m, const Size &size, const Values<double> &values)
//     {
//         m.values(size.get(0), size.get(1), values.value());
//     }

//     void BLASBackend::build(Matrix &m, const Size &size, const Zeros & /*values*/)
//     {
//         build(m, size, Values<double>(0));
//     }

//     void BLASBackend::build(Vector &m, const Size &size, const Zeros & /*values*/)
//     {
//         build(m, size, Values<double>(0));
//     }

//     void BLASBackend::build(Vector &m, const Size &size, const LocalZeros & /*values*/)
//     {
//         build(m, size, Values<double>(0));
//     }

//     void BLASBackend::build(Vector &v, const Size &size, const Values<double> &values)
//     {
//         v = Vector(size.get(0), values.value());
//     }

//     void BLASBackend::mat_diag_shift(Matrix &left, const BLASBackend::Scalar diag_factor)
//     {
//         using std::min;
//         const SizeType n = min(left.rows(), left.cols());
//         for(SizeType i = 0; i != n; ++i) {
//             left.set(i, i, left.get(i, i) + diag_factor);
//         }
//     }

//     void BLASBackend::diag(Matrix &left, const Vector &right)
//     {
//         const SizeType n = right.size();
//         left.resize(n, n);
//         std::fill(left.entries().begin(), left.entries().end(), 0);

//         for(SizeType i = 0; i != n; ++i) {
//             left.set(i, i,right[i]);
//         }
//     }

//     void BLASBackend::diag(Vector &left, const CRSMatrix<Scalar>  &right)
//     {
//         const SizeType n = std::min(right.rows(), right.cols());
//         left.resize(n);

//         for (SizeType r = 0; r < n; ++r) {
//             const SizeType begin = right.rowptr()[r];
//             const SizeType end   = right.rowptr()[r + 1];

//             left[r] = 0.;

//             for (SizeType ind = begin; ind != end; ++ind) {
//                 const SizeType c = right.colindex()[ind];
//                 if(c == r) {
//                     left[r] = right.entries()[ind];
//                     break;
//                 }
//             }
//         }
//     }

//     void BLASBackend::diag(Vector &left, const Matrix &right)
//     {
//         const SizeType n = std::min(right.rows(), right.cols());
//         left.resize(n);

//         for(SizeType i = 0; i != n; ++i) {
//             left[i] = right.get(i, i);
//         }
//     }

//     void BLASBackend::diag(Matrix &left, const Matrix &right)
//     {
//         const SizeType n = std::min(right.rows(), right.cols());
//         left.resize(n, n);

//         std::fill(left.entries().begin(), left.entries().end(), 0);

//         for(SizeType i = 0; i != n; ++i) {
//             left.set(i, i, right.get(i, i));
//         }
//     }

//     //-------------------------------------------------------------------------------------- //

//     void BLASBackend::apply_binary(Vector &result, const Reciprocal<Scalar> &reciprocal, const Vector &vec)
//     {
//         result.resize(vec.size());
//         for(SizeType i = 0; i < vec.size(); ++i) {
//             result[i] = reciprocal.numerator() / vec[i];
//         }
//     }

//     void BLASBackend::kronecker_product(Matrix &result, const Vector &left, const Vector &right)
//     {
//         result.resize(left.size(), right.size());
//         for (SizeType i = 0; i < left.size(); ++i) {
//             for (SizeType j = 0; j < right.size(); ++j) {
//                 result.set(i, j, left[i] * right[j]);
//             }
//         }
//     }

//     void BLASBackend::assign_transposed(Matrix &out, const Matrix &in)
//     {
//         out.resize(in.cols(), in.rows());
//         for(SizeType i = 0; i < in.rows(); ++i) {
//             for(SizeType j = 0; j < in.cols(); ++j) {
//                 out.set(j, i, in.get(i, j));
//             }
//         }
//     }

//     BLASBackend::Scalar BLASBackend::trace(const Matrix &in)
//     {
//         const SizeType n = std::min(in.rows(), in.cols());

//         BLASBackend::Scalar ret = 0;
//         for(SizeType i = 0; i < n; ++i) {
//             ret += in.get(i, i);
//         }

//         return ret;
//     }

//     void BLASBackend::apply_binary(Vector &result, const Vector &vec, const Multiplies &, const Scalar value)
//     {
//         axpy(result, value, vec);
//     }

//     void BLASBackend::apply_binary(Vector &result, const Vector &vec, const Divides &, const Scalar value)
//     {
//         axpy(result, 1./value, vec);
//     }

//     void BLASBackend::apply_binary(Vector &result, const Vector &left, const Minus &, const Vector &right)
//     {
//         result = left;
//         axpy(result, -1., right);
//     }

//     void BLASBackend::apply_binary(Vector &result, const Vector &left, const Plus &, const Vector &right)
//     {
//         result = left;
//         axpy(result, 1., right);
//     }

//     void BLASBackend::apply_binary(Vector &result, const Matrix &left, const Multiplies &, const Vector &right)
//     {
//         gemv(result, 0., 1., false, left, right);
//     }

//     void BLASBackend::apply_binary(Vector &result, const CRSMatrix<Scalar> &left, const Multiplies &, const Vector &right)
//     {
//         gemv(result, 0., 1., false, left, right);
//     }

//     void BLASBackend::apply_binary(Matrix &result, const Matrix &left, const Plus &, const Matrix &right)
//     {
//         result = left;
//         axpy(result, 1., right);
//     }

//     void BLASBackend::apply_binary(Matrix &result, const Matrix &left, const Minus & , const Matrix &right)
//     {
//         result = left;
//         axpy(result, -1., right);
//     }

//     void BLASBackend::apply_binary(Matrix &result, const Matrix &mat, const Multiplies &, const Scalar &value)
//     {
//         axpy(result, value, mat);
//     }

//     void BLASBackend::apply_binary(Matrix &result, const Matrix &mat, const Divides &, const Scalar &value)
//     {
//         axpy(result, 1./value, mat);
//     }

//     void BLASBackend::apply_binary(Matrix &result, const Matrix &left, const Multiplies &, const Matrix &right)
//     {
//         gemm(result, 0., 1, false, left, false, right);
//     }


//     bool BLASBackend::read(const std::string &path, CRSMatrix<Scalar> &mat)
//     {
//         SizeType offset = -1;
//         std::ifstream is(path.c_str());

//         if(!is.good()) {
//             is.close();
//             return false;
//         }

//         mat.clear();

//         std::string line;
//         std::getline(is, line);

//         SizeType rows, nnz;
//         int nvals = std::sscanf(line.c_str(), "%ld %ld", &rows, &nnz);

//         if(nvals != 2) {
//             std::cerr << "bad format for line 0: " << line << std::endl;
//         }

//         mat.initialize(rows, rows, nnz);

//         mat.assembly_begin();

//         SizeType line_num = 1;
//         while(is.good()) {
//             std::getline(is, line);

//             SizeType i, j;
//             Scalar val;

//             nvals = std::sscanf(line.c_str(), "%ld %ld %lg", &i, &j, &val);

//             if(nvals != 3) {
//                 std::cerr << "[Warning] bad format for line" << line_num << ": " << line << std::endl;
//                 continue;
//             }

//             mat.set(i + offset, j + offset, val);

//             if(line_num == nnz) break;
//             line_num++;
//         }

//         if(line_num != nnz) {
//             std::cerr << "[Warning] bad format, number of entries not equal to nnz" << std::endl;
//         }

//         mat.assembly_end();

//         is.close();
//         return false;
//     }

//     bool BLASBackend::read(const std::string &path, Vector &vec)
//     {
//         std::ifstream is(path.c_str());
//         vec.clear();

//         std::string line;
//         while(is.good()) {
//             std::getline(is, line);
//             vec.push_back(atof(line.c_str()));
//         }

//         is.close();
//         return false;
//     }

//     void BLASBackend::disp(const Vector &vec)
//     {
//         for(auto v : vec) {
//             std::cout << v << " ";
//         }

//         std::cout << std::endl;
//     }

//     void BLASBackend::disp(const Matrix &mat)
//     {
//         for(SizeType i = 0; i < mat.rows(); ++i) {
//             for(SizeType j = 0; j < mat.cols(); ++j) {
//                 std::cout << mat.get(i, j) << "\t";
//             }

//             std::cout << std::endl;
//         }
//     }

//     // void BLASBackend::disp(const CRSMatrix<Scalar> &mat)
//     // {
//     // 	disp(mat, std::cout);
//     // }

//     void BLASBackend::build_from_structure(Vector &lhs, const Vector &rhs)
//     {
//         lhs.resize(rhs.size());
//         std::fill(lhs.begin(), lhs.end(), 0.);
//     }

// }

// #endif //utopia_utopia_BLASBACKEND_CPP


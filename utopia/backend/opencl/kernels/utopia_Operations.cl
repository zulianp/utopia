#ifndef UTOPIA_OPERATIONS_CL
#define UTOPIA_OPERATIONS_CL

//TODO
//matrix_vector_multiplication_at_row(index, rows, columns, arg_1, arg_2)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#ifndef Scalar
typedef double Scalar;
#endif

#ifndef SizeType
typedef int SizeType;
#endif


Scalar pow2(Scalar macro_x);
Scalar pow2(Scalar macro_x) { return macro_x * macro_x; }

#define plus(macro_left, macro_right) (macro_left + macro_right)
#define minus(macro_left, macro_right) (macro_left - macro_right)
#define multiplies(macro_left, macro_right) (macro_left * macro_right)
#define divides(macro_left, macro_right) (macro_left / macro_right)
#define generic_copy(n, src, dest) { for(SizeType macro_index = 0; macro_index < (n); ++macro_index) { (dest)[macro_index] = (src)[macro_index]; } }

Scalar matrix_vector_multiplication_at_row(const SizeType r, const SizeType rows, const SizeType cols, __global const Scalar *mat, __global const Scalar *vec);

Scalar matrix_vector_multiplication_at_row(const SizeType r, const SizeType rows, const SizeType cols, __global const Scalar *mat, __global const Scalar *vec)
{

    const SizeType offset = r * cols;
    Scalar ret = 0;
    for(SizeType j = 0; j < cols; ++j) {
        ret += mat[offset + j] * vec[j];
    }

    return ret;
}

Scalar matrix_matrix_multiplication_at_entry(const SizeType i, const SizeType j, const SizeType rows, const SizeType cols, const SizeType right_cols, __global const Scalar *left, __global const Scalar *right);

Scalar matrix_matrix_multiplication_at_entry(const SizeType i, const SizeType j, const SizeType rows, const SizeType cols, const SizeType right_cols, __global const Scalar *left, __global const Scalar *right)
{
    const SizeType offset_left  = i * cols;
    // const SizeType offset_right = j * right_cols;

    Scalar result = 0;
    for(SizeType k = 0; k < cols; ++k) {
        result += left[offset_left + k] * right[j + k * right_cols];
    }

    return result;
}

int get_index(const int r, const int c, const int rows, const int columns);
int get_index_transposed(const int r, const int c, const int rows, const int columns);

int get_index(const int r, const int c, const int rows, const int columns)
{
    return r * columns + c;
}

int get_index_transposed(const int r, const int c, const int rows, const int columns)
{
    return c * rows + r;
}


#endif //UTOPIA_OPERATIONS_CL
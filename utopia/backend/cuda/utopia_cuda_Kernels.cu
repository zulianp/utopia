//#ifdef WITH_CUDA
//
//#ifndef UTOPIA_CUDA_KERNELS_H
//#define UTOPIA_CUDA_KERNELS_H
//
//namespace utopia {
//    namespace cuda {
//        template<typename T, class SizeType, class Operation>
//        __global__ void apply(const T *left, const T *right, T *result, SizeType n) {
//            Operation op;
//            const unsigned int nThreads = gridDim.x * blockDim.x;
//            for (unsigned int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += nThreads) {
//                result[i] = op(left[i], right[i]);
//            }
//        }
//
//        template<typename T>
//        class Add {
//        public:
//            __host__ __device__ T operator()(const T &left, const T &right) const {
//                return left + right;
//            }
//        };
//
//        template<typename T, class SizeType>
//        __host__ void sum(const T *left, const T *right, T *result, SizeType n) {
//            apply<T, SizeType, Add<T> > <<<n, 1>>>(left, right, n);
//        }
//    }
//}
//
//#endif //UTOPIA_CUDA_KERNELS_H
//#endif //WITH_CUDA
//
//
//
//

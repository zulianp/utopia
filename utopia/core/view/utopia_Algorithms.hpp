#ifndef UTOPIA_ALGORITHMS_HPP
#define UTOPIA_ALGORITHMS_HPP

#include "utopia_Base.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

#ifdef WITH_TRILINOS
#include <Kokkos_ArithTraits.hpp>
#endif

namespace utopia {
    namespace device {

#ifdef KOKKOS_INLINE_FUNCTION

        template<typename T>
        UTOPIA_INLINE_FUNCTION T abs(const T &v)
        {
            return Kokkos::Details::ArithTraits<T>::abs(v);
        }

        template<typename T>
        UTOPIA_INLINE_FUNCTION T min(const T &left, const T &right)
        {
            return left < right ? left : right;
        }

        template<typename T>
        UTOPIA_INLINE_FUNCTION T max(const T &left, const T &right)
        {
            return left > right ? left : right;
        }

        template<typename T>
        UTOPIA_INLINE_FUNCTION T sqrt(const T &value) {
            return Kokkos::Details::ArithTraits<T>::sqrt(value);
        }
    
#else

        template<typename T>
        inline T abs(const T &v)
        {
            return std::abs(v);
        }

        template<typename T>
        inline T min(const T &left, const T &right)
        {
            return std::min(left, right);
        }

        template<typename T>
        inline T max(const T &left, const T &right)
        {
            return std::max(left, right);
        }

        inline Scalar sqrt(const Scalar &value) {
            return std::sqrt(value);
        }
    
#endif //KOKKOS_INLINE_FUNCTION


        //FIXME use impl and specialize based on Input type
        template<class InputIterator, class OutputIterator>
        UTOPIA_INLINE_FUNCTION void copy(
            const InputIterator &in_begin,
            const InputIterator &in_end,
            const OutputIterator &out_begin)
        {
            auto in_it = in_begin;
            auto out_it = out_begin;
            
            for(; in_it != in_end; ++in_it, ++out_it)
            {
                *out_it = *in_it;
            }
        }

        template<class InputArray, class OutputArray>
        class Copy {
        public:

            UTOPIA_INLINE_FUNCTION static void apply(const InputArray &in, OutputArray &out)
            {
                UTOPIA_DEVICE_ASSERT(in.size() == out.size());

                const SizeType n = in.size();
                for(SizeType i = 0; i < n; ++i) {
                    out[i] = in[i];
                }
            }

        };

        template<class InputArray, class OutputArray>
        UTOPIA_INLINE_FUNCTION void copy(const InputArray &in, OutputArray &out)
        {
            Copy<InputArray, OutputArray>::apply(in, out);
        }

        template<class Array>
        class Fill {
        public:

            template<typename Value>
            UTOPIA_INLINE_FUNCTION static void apply(const Value &val, Array &in_out)
            {
                const SizeType n = in_out.size();
                for(SizeType i = 0; i < n; ++i) {
                    in_out[i] = val;
                }
            }
        };

        template<typename Value, class Array>
        UTOPIA_INLINE_FUNCTION static void fill(const Value &val, Array &in_out)
        {
            Fill<Array>::apply(val, in_out);
        }

        template<class Array, typename Value>
        UTOPIA_INLINE_FUNCTION static void apply(const Value &val, Array &out)
        {
            Fill<Array>::apply(val, out);
        }

        template<class NDimArray>
        class Extent {
        public:
            template<typename S>
            UTOPIA_INLINE_FUNCTION static auto apply(const NDimArray &in, const S &s) -> decltype(in.extent(s))
            {
                return in.extent(s);
            }  
        };

        template<class NDimArray, typename S>
        UTOPIA_INLINE_FUNCTION auto extent(const NDimArray &in, const S &s) -> decltype( Extent<NDimArray>::apply(in, s) )
        {
            return Extent<NDimArray>::apply(in, s);
        }

        template<class Array>
        class Scale {
        public:
            template<typename Scalar>
            UTOPIA_INLINE_FUNCTION static void apply(const Scalar &alpha, Array &in_out)
            {
                const SizeType n = in_out.size();
                for(SizeType i = 0; i < n; ++i) {
                    in_out[i] *= alpha;
                }
            }  
        };

        template<class Array, typename Scalar>
        UTOPIA_INLINE_FUNCTION void scale(const Scalar &alpha, Array &in_out) 
        {
            Scale<Array>::apply(alpha, in_out);
        }

        template<class X, class Y>
        class AXPY {
        public:
            template<typename Scalar>
            UTOPIA_INLINE_FUNCTION static void apply(const Scalar &alpha, const X &x, Y &y)
            {
                UTOPIA_DEVICE_ASSERT(x.size() == y.size());

                const SizeType n = x.size();
                for(SizeType i = 0; i < n; ++i) {
                    y[i] += alpha * x[i];
                }
            }  
        };

        template<class X, class Y, class Scalar>
        UTOPIA_INLINE_FUNCTION void axpy(const Scalar &alpha, const X &x, Y &y)
        {
            AXPY<X, Y>::apply(alpha, x, y);
        }

        template<class X, class Y>
        class Dot {
        public:

            UTOPIA_INLINE_FUNCTION static auto apply(const X &x, const Y &y) -> decltype(x[0]*y[0])
            {
                UTOPIA_DEVICE_ASSERT(x.size() == y.size());
                const SizeType n = x.size();
                if(n == 0) return 0;

                auto ret = x[0]*y[0];

                for(SizeType i = 1; i < n; ++i) {
                    ret += x[i] * y[i];
                }

                return ret;
            }  
        };

        template<class X, class Y>
        UTOPIA_INLINE_FUNCTION auto dot(const X &x, const Y &y) -> decltype(Dot<X, Y>::apply(x, y))
        {
            return Dot<X, Y>::apply(x, y);
        }

        template<class Matrix, class InVector, class OutVector>
        class MV {
        public:

            UTOPIA_INLINE_FUNCTION static void apply(const Matrix &m, const InVector &in, OutVector &out)
            {
                UTOPIA_DEVICE_ASSERT(in.size()  == extent(m, 1));
                UTOPIA_DEVICE_ASSERT(out.size() == extent(m, 0));

                const SizeType r = out.size();
                const SizeType c = in.size();

                for(SizeType i = 0; i < r; ++i) {
                    out[i] = 0;
                    for(SizeType j = 0; j < c; ++j) {
                        out[i] += m(i, j) * in[j];
                    }
                }
            }  
        };

        template<class Matrix, class InVector, class OutVector>
        UTOPIA_INLINE_FUNCTION void mv(const Matrix &m, const InVector &in, OutVector &out)
        {
            MV<Matrix, InVector, OutVector>::apply(m, in, out);
        }

        template<class Matrix, class InMatrix, class OutMatrix>
        class MM {
        public:

            UTOPIA_INLINE_FUNCTION static void apply(const Matrix &m, const InMatrix &in, OutMatrix &out)
            {
                UTOPIA_DEVICE_ASSERT(in.size() == extent(m, 1));
                UTOPIA_DEVICE_ASSERT(out.size() == extent(m, 0));

                const SizeType r = extent(m, 0);
                const SizeType l = extent(m, 1);
                const SizeType c = extent(in, 1);

                out.set(0.0);

                for(SizeType i = 0; i < r; ++i) {
                    for(SizeType j = 0; j < c; ++j) {
                        for(SizeType k = 0; k < l; ++k) {
                            out(i, j) += m(i, k) * in(k, j);
                        }
                    }
                }
            }  
        };

        template<class Matrix, class InMatrix, class OutMatrix>
        UTOPIA_INLINE_FUNCTION void mm(const Matrix &m, const InMatrix &in, OutMatrix &out)
        {
            MM<Matrix, InMatrix, OutMatrix>::apply(m, in, out);
        }

        template<class X, class Y>
        class ApproxEq {
        public:
            template<typename Scalar>
            UTOPIA_INLINE_FUNCTION static bool apply(const X &x, const Y &y, const Scalar &tol)
            {
                UTOPIA_DEVICE_ASSERT(x.size() == y.size());

                const SizeType n = x.size();
                for(SizeType i = 0; i < n; ++i) {
                    if(device::abs(x[i] - y[i]) > tol) {
                        return false;
                    }
                }

                return true;
            }  
        };

        template<class X, class Y, class Scalar>
        UTOPIA_INLINE_FUNCTION bool approxeq(const X &x, const Y &y, const Scalar &tol)
        {
            return ApproxEq<X, Y>::apply(x, y, tol);
        }

        template<class Scalar>
        UTOPIA_INLINE_FUNCTION bool approxeq(const Scalar &x, const Scalar &y, const Scalar &tol)
        {
            return device::abs(x - y) <= tol;
        }

        template<class Array>
        class Shift {
        public:
            template<typename Scalar>
            UTOPIA_INLINE_FUNCTION static void apply(const Scalar &alpha, Array &in_out)
            {
                const SizeType n = in_out.size();
                for(SizeType i = 0; i < n; ++i) {
                    in_out[i] += alpha;
                }
            }  
        };

        template<class Array, typename Scalar>
        UTOPIA_INLINE_FUNCTION void shift(const Scalar &alpha, Array &in_out) 
        {
            Shift<Array>::apply(alpha, in_out);
        }
    }

}

#endif //UTOPIA_ALGORITHMS_HPP

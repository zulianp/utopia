#ifndef UTOPIA_SIMD_OPS_HPP
#define UTOPIA_SIMD_OPS_HPP

#include "utopia_Base.hpp"

namespace utopia {
    namespace simd_v2 {

        template <class Impl>
        class Ops {};

        template <typename T, int Size_>
        class Auto {
        public:
            static constexpr size_t MemoryAlignment = alignof(T);
            static constexpr int Size = Size_;
            T data[Size];

            constexpr Auto(const T* data) {
                for (int i = 0; i < Size; ++i) {
                    this->data[i] = data[i];
                }
            }

            constexpr Auto(const T& value = 0) {
                for (int i = 0; i < Size; ++i) {
                    this->data[i] = value;
                }
            }

            void copy(const T* in) {
                for (int i = 0; i < Size; ++i) {
                    this->data[i] = in[i];
                }
            }

            inline Auto& operator=(const T& val) {
                for (int i = 0; i < Size; ++i) {
                    this->data[i] = val;
                }
                return *this;
            }

            inline Auto& operator+=(const Auto& val) {
                for (int i = 0; i < Size; ++i) {
                    this->data[i] += val.data[i];
                }
                return *this;
            }

            inline T sum() const {
                auto val = data[0];
                for (int i = 1; i < Size; ++i) {
                    val += this->data[i];
                }
                return val;
            }
        };

        template <typename T, int Lanes>
        class Ops<Auto<T, Lanes>> {
        public:
            using SIMDType = simd_v2::Auto<T, Lanes>;

            UTOPIA_INLINE_FUNCTION SIMDType construct(const T* data) {
                SIMDType ret;
                ret.copy(data);
                return ret;
            }

            UTOPIA_INLINE_FUNCTION static void add(const T* l, const T* r, T* result) {
                for (int i = 0; i < Lanes; ++i) {
                    result[i] = l[i] + l[r];
                }
            }

            UTOPIA_INLINE_FUNCTION static void in_place_add(const T* val, T* result) {
                for (int i = 0; i < Lanes; ++i) {
                    result[i] += val[i];
                }
            }

            UTOPIA_INLINE_FUNCTION static void subtract(const T* l, const T* r, T* result) {
                for (int i = 0; i < Lanes; ++i) {
                    result[i] = l[i] - l[r];
                }
            }

            UTOPIA_INLINE_FUNCTION static void in_place_subtract(const T* val, T* result) {
                for (int i = 0; i < Lanes; ++i) {
                    result[i] -= val[i];
                }
            }

            UTOPIA_INLINE_FUNCTION static void divide(const T* left, const T* right, T* result) {
                for (int i = 0; i < Lanes; ++i) {
                    result[i] = left[i] / right[i];
                }
            }

            UTOPIA_INLINE_FUNCTION static void in_place_divide(const T* denominator, T* result) {
                divide(result, denominator, result);
            }

            UTOPIA_INLINE_FUNCTION static void scale(const T& val, T* result) {
                for (int l = 0; l < Lanes; ++l) {
                    result[l] *= val;
                }
            }

            UTOPIA_INLINE_FUNCTION static void scale(const SIMDType& val, T* result) {
                for (int l = 0; l < Lanes; ++l) {
                    result[l] *= val[l];
                }
            }

            template <int N>
            UTOPIA_INLINE_FUNCTION static void dot(const T* left, const T* right, SIMDType& ret) {
                for (int l = 0; l < Lanes; ++l) {
                    ret.data[l] = left[l] * right[l];
                }

                for (int i = 1; i < N; ++i) {
                    for (int l = 0; l < Lanes; ++l) {
                        ret.data[l] += left[i * Lanes + l] * right[i * Lanes + l];
                    }
                }
            }

            UTOPIA_INLINE_FUNCTION void add_scale(const T& factor, const T* left, const T* right, T* result) {
                for (int l = 0; l < Lanes; ++l) {
                    result[l] = factor * (left[l] + right[l]);
                }
            }

            UTOPIA_INLINE_FUNCTION void add_scale_store(const T& factor, T* left, T* right) {
                for (int l = 0; l < Lanes; ++l) {
                    left[l] = factor * (left[l] + right[l]);
                    right[l] = left[l];
                }
            }
        };

    }  // namespace simd_v2

}  // namespace utopia

#ifdef UTOPIA_WITH_VC
#include <Vc/Vc>

namespace utopia {
    namespace simd_v2 {

        template <typename T>
        class Ops<Vc::Vector<T>> {
        public:
            using SIMDType = Vc::Vector<T>;
            static constexpr int Lanes = SIMDType::Size;

            inline static SIMDType construct(const T* v) { return load(v); }
            inline static SIMDType load(const T* v) { return SIMDType(v, Vc::Aligned); }
            inline static void store(const SIMDType& in, T* out) { in.store(out, Vc::Aligned); }
            inline static void add(const T* l, const T* r, T* result) { store(load(l) + load(r), result); }
            inline static void in_place_add(const T* val, T* result) { add(result, val, result); }
            inline static void subtract(const T* l, const T* r, T* result) { store(load(l) - load(r), result); }
            inline static void in_place_subtract(const T* val, T* result) { subtract(result, val, result); }

            inline static void divide(const T* left, const T* right, T* result) {
                store(load(left) / load(right), result);
            }

            inline static void divide(const T* left, const T& right, T* result) { store(load(left) / right, result); }

            inline static void in_place_divide(const T* denominator, T* result) { divide(result, denominator, result); }
            inline static void in_place_divide(const T& denominator, T* result) { divide(result, denominator, result); }

            inline static void scale(const T& val, T* result) { store(load(result) * val, result); }
            inline static void scale(const SIMDType& val, T* result) { store(load(result) * val, result); }

            template <int N>
            inline static void dot(const T* l, const T* r, SIMDType& ret) {
                ret = load(&l[0]) * load(&r[0]);

                for (int i = 1; i < N; ++i) {
                    ret += load(&l[i * Lanes]) * load(&r[i * Lanes]);
                }
            }

            inline void add_scale(const T& factor, const T* l, const T* r, T* result) {
                store(factor * (load(l) + load(r)), result);
            }

            inline void add_scale_store(const T& factor, T* l, T* r) {
                auto res = factor * (load(l) + load(r));
                store(res, l);
                store(res, r);
            }
        };

    }  // namespace simd_v2
}  // namespace utopia

#endif  // UTOPIA_WITH_VC

#endif  // UTOPIA_SIMD_OPS_HPP

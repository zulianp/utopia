#ifndef UTOPIA_SIMD_VECTOR_HPP
#define UTOPIA_SIMD_VECTOR_HPP

#include "utopia_simd_Ops.hpp"

namespace utopia {

    namespace simd_v2 {

        /**
         * Storage will be N * Lanes
         */
        template <typename T, int N_, typename SIMDType_ = Vc::Vector<T>>
        class
            // __attribute__((__packed__))
            Vector final {
        public:
            using SIMDType = SIMDType_;
            static constexpr int N = N_;
            static constexpr int Lanes = SIMDType::Size;
            static constexpr int Size = N * Lanes;
            using Ops = utopia::simd_v2::Ops<SIMDType>;

            // enum { StoreAs = UTOPIA_BY_REFERENCE };

            // static constexpr size_t data_alignment() { return SIMDType::MemoryAlignment; }
            // alignas(data_alignment())
            T data[N * Lanes];

            template <typename Factor>
            inline void scale(const Factor &factor) {
                for (int i = 0; i < N; ++i) {
                    Ops::scale(factor, block(i));
                }
            }

            template <typename Factor>
            inline Vector operator*(const Factor &factor) {
                Vector ret = *this;
                ret.scale(factor);
                return ret;
            }

            inline SIMDType operator[](const int idx) const { return Ops::construct(block(idx)); }

            inline T &ref(const int component, const int lane) {
                assert(component < N);
                assert(lane < Lanes);
                return data[component * Lanes + lane];
            }

            inline constexpr const T &get(const int component, const int lane) const {
                assert(component < N);
                assert(lane < Lanes);
                return data[component * Lanes + lane];
            }

            inline constexpr const T *block(const int idx) const {
                assert(idx < N);
                assert(idx >= 0);

                return &data[idx * Lanes];
            }

            inline T *block(const int idx) {
                assert(idx < N);
                assert(idx >= 0);
                return &data[idx * Lanes];
            }

            inline constexpr Vector &operator*=(const T &other) {
                scale(other);
                return *this;
            }

            inline constexpr Vector &operator/=(const T &other) {
                scale(1. / other);
                return *this;
            }

            inline constexpr Vector &operator/=(const Vector &other) {
                for (int i = 0; i < N; ++i) {
                    Ops::in_place_divide(other.block(i), block(i));
                }

                return *this;
            }

            inline constexpr Vector &operator+=(const Vector &other) {
                for (int i = 0; i < N; ++i) {
                    Ops::in_place_add(other.block(i), block(i));
                }

                return *this;
            }

            inline constexpr Vector &set(const int idx, const SIMDType &val) {
                Ops::store(val, block(idx));
                return *this;
            }

            inline constexpr Vector &add(const int idx, const SIMDType &val) {
                Ops::in_place_add(val, block(idx));
                return *this;
            }
            inline constexpr Vector &add(const int idx, const T &val) {
                Ops::in_place_add(val, block(idx));
                return *this;
            }

            inline constexpr Vector &divide(const int idx, const T &val) {
                Ops::in_place_divide(val, block(idx));
                return *this;
            }

            inline constexpr Vector &operator-=(const Vector &other) {
                for (int i = 0; i < N; ++i) {
                    Ops::in_place_subtract(other.block(i), block(i));
                }

                return *this;
            }

            inline void axpy(const int idx, const T &alpha, const SIMDType &x) {
                auto alpha_x = alpha * x;
                add(idx, alpha_x);
            }

            inline friend constexpr Vector operator+(const Vector &l, const Vector &r) {
                Vector ret = l;
                ret += r;
                return ret;
            }

            inline void dot(const Vector &other, SIMDType &result) const {
                Ops::template dot<N>(data, other.data, result);
            }

            friend inline constexpr SIMDType dot(const Vector &l, const Vector &r) {
                SIMDType ret;
                l.dot(r, ret);
                return ret;
            }

            friend inline constexpr SIMDType inner(const Vector &l, const Vector &r) {
                SIMDType ret;
                l.dot(r, ret);
                return ret;
            }

            friend void disp(const Vector &v, std::ostream &os = std::cout) {
                for (int i = 0; i < Size; i += Lanes) {
                    os << "[ ";
                    for (int l = 0; l < Lanes; ++l) {
                        os << v.data[i + l] << " ";
                    }

                    os << "]";
                    os << "\n";
                }
            }

            void set(const T &val) {
                for (int i = 0; i < Size; ++i) {
                    data[i] = val;
                }
            }
        };
    }  // namespace simd_v2

    // template <typename Array, typename T2, int Cols>
    // inline auto operator*(const TensorView<Array, 2> &mat, const Vector<T2, Cols> &x)
    //     -> Vector<decltype(mat(0, 0) * T2()), Array::Rows> {
    //     using RetT = decltype(mat(0, 0) * T2());

    //     Vector<RetT, Array::Rows> ret;

    //     for (int i = 0; i < static_cast<int>(Array::Rows); ++i) {
    //         ret[i] = Zero<RetT>::value();

    //         for (int j = 0; j < Cols; ++j) {
    //             ret[i] += mat(i, j) * x[j];
    //         }
    //     }

    //     return ret;
    // }

    template <typename T, int N, typename SIMDType_>
    class Traits<simd_v2::Vector<T, N, SIMDType_>> {
    public:
        using Primitive = T;
        using Scalar = SIMDType_;
        using SizeType = int;
        static const int Order = 1;
    };

}  // namespace utopia

#endif  // UTOPIA_SIMD_VECTOR_HPP

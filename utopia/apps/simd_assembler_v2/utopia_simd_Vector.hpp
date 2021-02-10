#ifndef UTOPIA_SIMD_VECTOR_HPP
#define UTOPIA_SIMD_VECTOR_HPP

namespace utopia {

    namespace simd_v2 {

        /**
         * Storage will be N * Lanes
         */
        template <typename T, int N_, int Lanes_ = Vc::Vector<T>::Size>
        class Vector final {
        public:
            using SIMDType = Vc::Vector<T>;
            static constexpr int N = N_;
            static constexpr int Lanes = Lanes_;
            static constexpr int Size = N * Lanes;

            enum { StoreAs = UTOPIA_BY_REFERENCE };

            static constexpr size_t data_alignment() { return SIMDType::MemoryAlignment; }
            alignas(data_alignment()) T data[N * Lanes];

            template <typename Factor>
            inline void scale(const Factor &factor) {
                SIMDType temp;
                for (int i = 0; i < N; ++i) {
                    temp = get(i) * factor;
                    set(i, temp);
                }
            }

            template <typename Factor>
            inline Vector operator*(const Factor &factor) {
                Vector ret = *this;
                ret.scale(factor);
                return ret;
            }

            inline constexpr const T &operator()(const int idx) const { return data[idx]; }

            inline constexpr SIMDType get(const int idx) const { return SIMDType(&data[idx * Lanes], Vc::Aligned); }
            inline void set(const int idx, const SIMDType &v) { v.store(&data[idx * Lanes], Vc::Aligned); }

            // inline T &operator()(const int idx) { return data[idx]; }

            // inline constexpr const T &operator[](const int idx) const { return data[idx]; }
            // inline T &operator[](const int idx) { return data[idx]; }

            // inline constexpr Vector operator+=(const Vector &other) {
            //     x() += other.x();
            //     y() += other.y();
            //     z() += other.z();
            //     return *this;
            // }

            // friend inline constexpr T dot(const Vector &l, const Vector &r) {
            //     return l.x() * r.x() + l.y() * r.y() + l.z() * r.z();
            // }

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
                for (auto &d : data) d = val;
            }
        };
    }  // namespace simd_v2

    // template <typename T, int Dim>
    // inline T inner(const Vector<T, Dim> &l, const Vector<T, Dim> &r) {
    //     return dot(l, r);
    // }

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

    template <typename T, int N, int Lanes>
    class Traits<simd_v2::Vector<T, N, Lanes>> {
    public:
        using Scalar = T;
        using SIMDType = Vc::Vector<T>;
    };

}  // namespace utopia

#endif  // UTOPIA_SIMD_VECTOR_HPP

#ifndef UTOPIA_SIMD_QUADRATURE_HPP
#define UTOPIA_SIMD_QUADRATURE_HPP

#include "utopia_DeviceExpression.hpp"
#include "utopia_Views.hpp"

#include <Vc/Vc>
#include <vector>

namespace utopia {

    namespace simd {
        template <typename T, int Dim, typename...>
        class Vector {};
    }  // namespace simd

    template <typename T, int Dim, typename... Args>
    class Traits<simd::Vector<T, Dim, Args...>> {
    public:
        using Scalar = T;
    };

    namespace simd {
        template <typename T>
        T sum(const T v) {
            return v;
        }

        template <typename T>
        class Vector<T, 2> final /*: public DeviceExpression<Vector<T, 3>>*/ {
        public:
            T data_[2] = {Zero<T>::value(), Zero<T>::value()};

            inline T &x() { return data_[0]; }
            inline constexpr const T &x() const { return data_[0]; }

            inline T &y() { return data_[1]; }
            inline constexpr const T &y() const { return data_[1]; }

            inline constexpr Vector operator*(const T &scale) { return {x() * scale, y() * scale}; }

            inline constexpr const T &operator()(const int idx) const { return data_[idx]; }
            inline T &operator()(const int idx) { return data_[idx]; }

            inline constexpr const T &operator[](const int idx) const { return data_[idx]; }
            inline T &operator[](const int idx) { return data_[idx]; }

            inline constexpr Vector operator+=(const Vector &other) {
                x() += other.x();
                y() += other.y();

                return *this;
            }

            friend inline constexpr T dot(const Vector &l, const Vector &r) { return l.x() * r.x() + l.y() * r.y(); }

            friend void disp(const Vector &v, std::ostream &os = std::cout) { os << v.x() << " " << v.y() << "\n"; }
        };

        template <typename T>
        class Vector<T, 3> final /*: public DeviceExpression<Vector<T, 3>>*/ {
        public:
            T data_[3] = {Zero<T>::value(), Zero<T>::value(), Zero<T>::value()};

            inline T &x() { return data_[0]; }
            inline constexpr const T &x() const { return data_[0]; }

            inline T &y() { return data_[1]; }
            inline constexpr const T &y() const { return data_[1]; }

            inline T &z() { return data_[2]; }
            inline constexpr const T &z() const { return data_[2]; }

            inline constexpr Vector operator*(const T &scale) { return {x() * scale, y() * scale, z() * scale}; }

            inline constexpr const T &operator()(const int idx) const { return data_[idx]; }
            inline T &operator()(const int idx) { return data_[idx]; }

            inline constexpr const T &operator[](const int idx) const { return data_[idx]; }
            inline T &operator[](const int idx) { return data_[idx]; }

            inline constexpr Vector operator+=(const Vector &other) {
                x() += other.x();
                y() += other.y();
                z() += other.z();
                return *this;
            }

            friend inline constexpr T dot(const Vector &l, const Vector &r) {
                return l.x() * r.x() + l.y() * r.y() + l.z() * r.z();
            }

            friend void disp(const Vector &v, std::ostream &os = std::cout) {
                os << v.x() << " " << v.y() << " " << v.z() << "\n";
            }
        };

        template <typename T, int Dim>
        class Quadrature {
        public:
            using Point = simd::Vector<T, Dim>;
            using ViewDevice = Quadrature;
            using ViewHost = Quadrature;

            Vc::vector<Point> points;
            Vc::vector<T> weights;

            inline int n_points() const { return weights.size(); }
            inline const Quadrature &view_device() const { return *this; }
            inline const Quadrature &view_host() const { return *this; }

            void resize(int n) {
                points.resize(n);
                weights.resize(n);
            }

            inline const Point &point(const int i) const { return points[i]; }
            inline const T &weight(const int i) const { return weights[i]; }
        };

        template <typename T>
        class Gauss {
        public:
            class Hex {
            public:
                inline static void set(Quadrature<T, 3> &q, const int idx, const T &x, const T &y, const T &z) {
                    q.points[idx][0] = x;
                    q.points[idx][1] = y;
                    q.points[idx][2] = z;
                }

                static bool get(const int order, Quadrature<T, 3> &q) {
                    switch (order) {
                        case 0:
                        case 1: {
                            q.resize(6);

                            set(q, 0, 0.0, 0.5, 0.5);
                            set(q, 1, 0.5, 0.0, 0.5);
                            set(q, 2, 0.5, 0.5, 0.0);
                            set(q, 3, 0.5, 0.5, 1.0);
                            set(q, 4, 0.5, 1.0, 0.5);
                            set(q, 5, 1.0, 0.5, 0.5);

                            q.weights[0] = 0.16666666666666666666666666666667;
                            q.weights[1] = 0.16666666666666666666666666666667;
                            q.weights[2] = 0.16666666666666666666666666666667;
                            q.weights[3] = 0.16666666666666666666666666666667;
                            q.weights[4] = 0.16666666666666666666666666666667;
                            q.weights[5] = 0.16666666666666666666666666666667;

                            return true;
                        }

                        case 2:
                        case 3:
                        case 4: {
                            q.resize(27);

                            set(q, 0, 0.112701665379258, 0.112701665379258, 0.112701665379258);
                            set(q, 1, 0.500000000000000, 0.112701665379258, 0.112701665379258);
                            set(q, 2, 0.887298334620742, 0.112701665379258, 0.112701665379258);
                            set(q, 3, 0.112701665379258, 0.500000000000000, 0.112701665379258);
                            set(q, 4, 0.500000000000000, 0.500000000000000, 0.112701665379258);
                            set(q, 5, 0.887298334620742, 0.500000000000000, 0.112701665379258);
                            set(q, 6, 0.112701665379258, 0.887298334620742, 0.112701665379258);
                            set(q, 7, 0.500000000000000, 0.887298334620742, 0.112701665379258);
                            set(q, 8, 0.887298334620742, 0.887298334620742, 0.112701665379258);
                            set(q, 9, 0.112701665379258, 0.112701665379258, 0.500000000000000);
                            set(q, 10, 0.500000000000000, 0.112701665379258, 0.500000000000000);
                            set(q, 11, 0.887298334620742, 0.112701665379258, 0.500000000000000);
                            set(q, 12, 0.112701665379258, 0.500000000000000, 0.500000000000000);
                            set(q, 13, 0.500000000000000, 0.500000000000000, 0.500000000000000);
                            set(q, 14, 0.887298334620742, 0.500000000000000, 0.500000000000000);
                            set(q, 15, 0.112701665379258, 0.887298334620742, 0.500000000000000);
                            set(q, 16, 0.500000000000000, 0.887298334620742, 0.500000000000000);
                            set(q, 17, 0.887298334620742, 0.887298334620742, 0.500000000000000);
                            set(q, 18, 0.112701665379258, 0.112701665379258, 0.887298334620742);
                            set(q, 19, 0.500000000000000, 0.112701665379258, 0.887298334620742);
                            set(q, 20, 0.887298334620742, 0.112701665379258, 0.887298334620742);
                            set(q, 21, 0.112701665379258, 0.500000000000000, 0.887298334620742);
                            set(q, 22, 0.500000000000000, 0.500000000000000, 0.887298334620742);
                            set(q, 23, 0.887298334620742, 0.500000000000000, 0.887298334620742);
                            set(q, 24, 0.112701665379258, 0.887298334620742, 0.887298334620742);
                            set(q, 25, 0.500000000000000, 0.887298334620742, 0.887298334620742);
                            set(q, 26, 0.887298334620742, 0.887298334620742, 0.887298334620742);

                            q.weights = {0.021433470507545, 0.034293552812071, 0.021433470507545, 0.034293552812071,
                                         0.054869684499314, 0.034293552812071, 0.021433470507545, 0.034293552812071,
                                         0.021433470507545, 0.034293552812071, 0.054869684499314, 0.034293552812071,
                                         0.054869684499314, 0.087791495198903, 0.054869684499314, 0.034293552812071,
                                         0.054869684499314, 0.034293552812071, 0.021433470507545, 0.034293552812071,
                                         0.021433470507545, 0.034293552812071, 0.054869684499314, 0.034293552812071,
                                         0.021433470507545, 0.034293552812071, 0.021433470507545};

                            return true;
                        }

                        default:
                            return false;
                    }
                }
            };
        };  // namespace simd

        template <typename T>
        class Gauss<Vc::Vector<T>> {
        public:
            using VectorType = Vc::Vector<T>;
            static const int block_size = VectorType::Size;

            class Hex {
            public:
                static bool get(const int order, simd::Quadrature<VectorType, 3> &q) {
                    Quadrature<T, 3> q_scalar;

                    if (!Gauss<T>::Hex::get(order, q_scalar)) {
                        return false;
                    }

                    int n_points = q_scalar.weights.size();
                    int n_blocks = n_points / block_size;
                    int n_defect = n_points - n_blocks * block_size;

                    int n_simd_blocks = n_blocks + (n_defect != 0);

                    q.resize(n_simd_blocks);

                    for (int i = 0; i < n_blocks; ++i) {
                        int q_offset = i * block_size;

                        for (int k = 0; k < block_size; ++k) {
                            int q_idx = q_offset + k;

                            q.points[i][0][k] = q_scalar.points[q_idx][0];
                            q.points[i][1][k] = q_scalar.points[q_idx][1];
                            q.points[i][2][k] = q_scalar.points[q_idx][2];
                            q.weights[i][k] = q_scalar.weights[q_idx];
                        }
                    }

                    if (n_defect) {
                        for (int k = 0; k < n_defect; ++k) {
                            int q_idx = n_blocks * block_size + k;

                            q.points[n_blocks][0][k] = q_scalar.points[q_idx][0];
                            q.points[n_blocks][1][k] = q_scalar.points[q_idx][1];
                            q.points[n_blocks][2][k] = q_scalar.points[q_idx][2];
                            q.weights[n_blocks][k] = q_scalar.weights[q_idx];
                        }

                        for (int k = n_defect; k < block_size; ++k) {
                            q.points[n_blocks][0][k] = Zero<T>::value();
                            q.points[n_blocks][1][k] = Zero<T>::value();
                            q.points[n_blocks][2][k] = Zero<T>::value();
                            q.weights[n_blocks][k] = Zero<T>::value();
                        }
                    }

                    return true;
                }
            };
        };

    }  // namespace simd
}  // namespace utopia
#endif  // UTOPIA_SIMD_QUADRATURE_HPP

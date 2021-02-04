#ifndef UTOPIA_SIMD_QUADRATURE_HPP
#define UTOPIA_SIMD_QUADRATURE_HPP

#include "utopia_DeviceExpression.hpp"
#include "utopia_Views.hpp"

#include "utopia_MultiVariateElement.hpp"
#include "utopia_UniformHex8.hpp"
#include "utopia_UniformQuad4.hpp"

#include "utopia_simd_common.hpp"

#include <Vc/Vc>
#include <vector>

namespace utopia {

    namespace simd {
        template <typename T, int Dim, typename...>
        class Vector {};

        template <typename T, int Rows, int Cols>
        class Matrix {
        public:
            static const int Size = Rows * Cols;
            T data_[Size];

            void set(const T &val) {
                for (int i = 0; i < Size; ++i) {
                    data_[i] = val;
                }
            }

            T &operator()(const int i, const int j) { return data_[i * Rows + j]; }
            const T &operator()(const int i, const int j) const { return data_[i * Rows + j]; }
        };
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
        inline T integrate(const T &v) {
            return sum(v);
        }

        template <typename T>
        inline T integrate(const Vc::Vector<T> &v) {
            return v.sum();
        }

        template <typename T>
        class Vector<T, 2> final /*: public DeviceExpression<Vector<T, 3>>*/ {
        public:
            T data_[2] = {simd::Zero<T>::value(), simd::Zero<T>::value()};

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
            T data_[3] = {simd::Zero<T>::value(), simd::Zero<T>::value(), simd::Zero<T>::value()};

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
        inline T inner(const Vector<T, Dim> &l, const Vector<T, Dim> &r) {
            return dot(l, r);
        }

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
            class Quad {
            public:
                inline static void set(Quadrature<T, 2> &q, const int idx, const T &x, const T &y) {
                    q.points[idx][0] = x;
                    q.points[idx][1] = y;
                }

                static bool get(const int order, Quadrature<T, 2> &q) {
                    switch (order) {
                        case 0:
                        case 1:
                        case 2: {
                            set(q, 0, 0.5, 0.5);
                            set(q, 1, 0.9304589153964795245728880523899, 0.5);
                            set(q, 2, 0.72780186391809642112479237299488, 0.074042673347699754349082179816666);
                            set(q, 3, 0.72780186391809642112479237299488, 0.92595732665230024565091782018333);
                            set(q, 4, 0.13418502421343273531598225407969, 0.18454360551162298687829339850317);
                            set(q, 5, 0.13418502421343273531598225407969, 0.81545639448837701312170660149683);

                            q.weights[0] = 0.28571428571428571428571428571428;
                            q.weights[1] = 0.10989010989010989010989010989011;
                            q.weights[2] = 0.14151805175188302631601261486295;
                            q.weights[3] = 0.14151805175188302631601261486295;
                            q.weights[4] = 0.16067975044591917148618518733485;
                            q.weights[5] = 0.16067975044591917148618518733485;
                            return true;
                        }
                        default:
                            return false;
                    }
                }
            };

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

            template <int Dim>
            static void vectorize(simd::Quadrature<T, Dim> &q_scalar, simd::Quadrature<VectorType, Dim> &q) {
                int n_points = q_scalar.weights.size();
                int n_blocks = n_points / block_size;
                int n_defect = n_points - n_blocks * block_size;

                int n_simd_blocks = n_blocks + (n_defect != 0);

                q.resize(n_simd_blocks);

                for (int i = 0; i < n_blocks; ++i) {
                    int q_offset = i * block_size;

                    for (int k = 0; k < block_size; ++k) {
                        int q_idx = q_offset + k;

                        for (int d = 0; d < Dim; ++d) {
                            q.points[i][d][k] = q_scalar.points[q_idx][d];
                        }

                        q.weights[i][k] = q_scalar.weights[q_idx];
                    }
                }

                if (n_defect) {
                    for (int k = 0; k < n_defect; ++k) {
                        int q_idx = n_blocks * block_size + k;

                        for (int d = 0; d < Dim; ++d) {
                            q.points[n_blocks][d][k] = q_scalar.points[q_idx][d];
                        }

                        q.weights[n_blocks][k] = q_scalar.weights[q_idx];
                    }

                    for (int k = n_defect; k < block_size; ++k) {
                        for (int d = 0; d < Dim; ++d) {
                            q.points[n_blocks][d][k] = Zero<T>::value();
                        }

                        q.weights[n_blocks][k] = Zero<T>::value();
                    }
                }
            }

            class Quad {
            public:
                static bool get(const int order, simd::Quadrature<VectorType, 2> &q) {
                    Quadrature<T, 2> q_scalar;

                    if (!Gauss<T>::Quad::get(order, q_scalar)) {
                        return false;
                    }

                    vectorize(q_scalar, q);
                    return true;
                }
            };

            class Hex {
            public:
                static bool get(const int order, simd::Quadrature<VectorType, 3> &q) {
                    Quadrature<T, 3> q_scalar;

                    if (!Gauss<T>::Hex::get(order, q_scalar)) {
                        return false;
                    }

                    vectorize(q_scalar, q);
                    return true;
                }
            };
        };

        template <typename...>
        class QuadratureDB {};

        template <typename T, typename QT>
        class QuadratureDB<UniformHex8<T>, QT> {
        public:
            template <class Q>
            static bool get(const int order, Q &q) {
                return Gauss<QT>::Hex::get(order, q);
            }
        };

        template <typename T, typename QT>
        class QuadratureDB<UniformQuad4<T>, QT> {
        public:
            template <class Q>
            static bool get(const int order, Q &q) {
                return Gauss<QT>::Quad::get(order, q);
            }
        };

        template <typename Elem, int NVar, typename QT>
        class QuadratureDB<MultiVariateElem<Elem, NVar>, QT> : public QuadratureDB<Elem, QT> {};

    }  // namespace simd
}  // namespace utopia
#endif  // UTOPIA_SIMD_QUADRATURE_HPP

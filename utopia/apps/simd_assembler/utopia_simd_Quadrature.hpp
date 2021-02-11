#ifndef UTOPIA_SIMD_QUADRATURE_HPP
#define UTOPIA_SIMD_QUADRATURE_HPP

#include "utopia_MultiVariateElement.hpp"
// #include "utopia_Tet4.hpp"
#include "utopia_Tri3.hpp"
#include "utopia_UniformHex8.hpp"
#include "utopia_UniformQuad4.hpp"

#include "utopia_simd_common.hpp"

#include <Vc/Vc>
#include <vector>

namespace utopia {

    namespace simd {

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
            inline static void set(Quadrature<T, 3> &q, const int idx, const T &x, const T &y, const T &z) {
                q.points[idx][0] = x;
                q.points[idx][1] = y;
                q.points[idx][2] = z;
            }

            inline static void set(Quadrature<T, 2> &q, const int idx, const T &x, const T &y) {
                q.points[idx][0] = x;
                q.points[idx][1] = y;
            }

            class Tri {
            public:
                static bool get(const int order, Quadrature<T, 2> &q) {
                    switch (order) {
                        case 0:
                        case 1:
                        case 2: {
                            q.resize(6);

                            set(q, 0, 0.5, 0.5);
                            set(q, 1, 0.5, 0.0);
                            set(q, 2, 0.0, 0.5);
                            set(q, 3, 1.0 / 6.0, 1.0 / 6.0);
                            set(q, 4, 1.0 / 6.0, 2.0 / 3.0);
                            set(q, 5, 2.0 / 3.0, 1.0 / 6.0);

                            q.weights[0] = 1.0 / 30.0;
                            q.weights[1] = 1.0 / 30.0;
                            q.weights[2] = 1.0 / 30.0;
                            q.weights[3] = 0.3;
                            q.weights[4] = 0.3;
                            q.weights[5] = 0.3;
                            return true;
                        }
                            // case 3:
                            // case 4:
                            // case 5:
                            // case 6: {
                            //     q.points = {{0.063089014491502228340331602870819,
                            //     0.063089014491502228340331602870819},
                            //                 {0.063089014491502228340331602870819,
                            //                 0.87382197101699554331933679425836}, {0.87382197101699554331933679425836,
                            //                 0.063089014491502228340331602870819},
                            //                 {0.24928674517091042129163855310702, 0.24928674517091042129163855310702},
                            //                 {0.24928674517091042129163855310702, 0.50142650965817915741672289378596},
                            //                 {0.50142650965817915741672289378596, 0.24928674517091042129163855310702},
                            //                 {0.053145049844816947353249671631398,
                            //                 0.31035245103378440541660773395655},
                            //                 {0.053145049844816947353249671631398,
                            //                 0.63650249912139864723014259441205}, {0.31035245103378440541660773395655,
                            //                 0.053145049844816947353249671631398},
                            //                 {0.31035245103378440541660773395655, 0.63650249912139864723014259441205},
                            //                 {0.63650249912139864723014259441205,
                            //                 0.053145049844816947353249671631398},
                            //                 {0.63650249912139864723014259441205,
                            //                 0.31035245103378440541660773395655}};

                            //     q.weights = {0.050844906370206816920936809106869,
                            //                  0.050844906370206816920936809106869,
                            //                  0.050844906370206816920936809106869,
                            //                  0.11678627572637936602528961138558,
                            //                  0.11678627572637936602528961138558,
                            //                  0.11678627572637936602528961138558,
                            //                  0.082851075618373575193553456420442,
                            //                  0.082851075618373575193553456420442,
                            //                  0.082851075618373575193553456420442,
                            //                  0.082851075618373575193553456420442,
                            //                  0.082851075618373575193553456420442,
                            //                  0.082851075618373575193553456420442};

                            //     return true;
                            // }
                            // case 7:
                            // case 8:
                            // case 9: {
                            //     q.points = {{0.333333333333333333333333333333333,
                            //     0.333333333333333333333333333333333},
                            //                 {0.48968251919873762778370692483619, 0.48968251919873762778370692483619},
                            //                 {0.48968251919873762778370692483619, 0.02063496160252474443258615032762},
                            //                 {0.02063496160252474443258615032762, 0.48968251919873762778370692483619},
                            //                 {0.43708959149293663726993036443535, 0.43708959149293663726993036443535},
                            //                 {0.43708959149293663726993036443535, 0.12582081701412672546013927112929},
                            //                 {0.12582081701412672546013927112929, 0.43708959149293663726993036443535},
                            //                 {0.18820353561903273024096128046733, 0.18820353561903273024096128046733},
                            //                 {0.18820353561903273024096128046733, 0.62359292876193453951807743906533},
                            //                 {0.62359292876193453951807743906533, 0.18820353561903273024096128046733},
                            //                 {0.044729513394452709865106589966276,
                            //                 0.044729513394452709865106589966276},
                            //                 {0.044729513394452709865106589966276,
                            //                 0.91054097321109458026978682006745}, {0.91054097321109458026978682006745,
                            //                 0.044729513394452709865106589966276},
                            //                 {0.74119859878449802069007987352342,
                            //                 0.036838412054736283634817598783385},
                            //                 {0.74119859878449802069007987352342, 0.22196298916076569567510252769319},
                            //                 {0.036838412054736283634817598783385,
                            //                 0.74119859878449802069007987352342},
                            //                 {0.036838412054736283634817598783385,
                            //                 0.22196298916076569567510252769319}, {0.22196298916076569567510252769319,
                            //                 0.74119859878449802069007987352342}, {0.22196298916076569567510252769319,
                            //                 0.036838412054736283634817598783385}};

                            //     q.weights = {0.097135796282798833819241982507289,
                            //                  0.031334700227139070536854831287209,
                            //                  0.031334700227139070536854831287209,
                            //                  0.031334700227139070536854831287209,
                            //                  0.077827541004774279316739356299404,
                            //                  0.077827541004774279316739356299404,
                            //                  0.077827541004774279316739356299404,
                            //                  0.079647738927210253032891774264045,
                            //                  0.079647738927210253032891774264045,
                            //                  0.079647738927210253032891774264045,
                            //                  0.025577675658698031261678798559000,
                            //                  0.025577675658698031261678798559000,
                            //                  0.025577675658698031261678798559000,
                            //                  0.043283539377289377289377289377289,
                            //                  0.043283539377289377289377289377289,
                            //                  0.043283539377289377289377289377289,
                            //                  0.043283539377289377289377289377289,
                            //                  0.043283539377289377289377289377289,
                            //                  0.043283539377289377289377289377289};

                            //     return true;
                            // }

                        default: {
                            assert(false);
                            return false;
                        }
                    }

                    return false;
                }
            };

            class Quad {
            public:
                static bool get(const int order, Quadrature<T, 2> &q) {
                    switch (order) {
                        case 0:
                        case 1:
                        case 2: {
                            q.resize(6);
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
                        default: {
                            assert(false);
                            return false;
                        }
                    }
                }
            };

            class Hex {
            public:
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

                        default: {
                            assert(false);
                            return false;
                        }
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

            class Tri {
            public:
                static bool get(const int order, simd::Quadrature<VectorType, 2> &q) {
                    Quadrature<T, 2> q_scalar;

                    if (!Gauss<T>::Tri::get(order, q_scalar)) {
                        return false;
                    }

                    vectorize(q_scalar, q);
                    return true;
                }
            };

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

        template <typename T, typename QT>
        class QuadratureDB<Tri3<T>, QT> {
        public:
            template <class Q>
            static bool get(const int order, Q &q) {
                return Gauss<QT>::Tri::get(order, q);
            }
        };

        // template <typename T, typename QT>
        // class QuadratureDB<Tet4<T>, QT> {
        // public:
        //     template <class Q>
        //     static bool get(const int order, Q &q) {
        //         return Gauss<QT>::Tet::get(order, q);
        //     }
        // };

        template <typename Elem, int NVar, typename QT>
        class QuadratureDB<MultiVariateElem<Elem, NVar>, QT> : public QuadratureDB<Elem, QT> {};

    }  // namespace simd
}  // namespace utopia

#endif  // UTOPIA_SIMD_QUADRATURE_HPP

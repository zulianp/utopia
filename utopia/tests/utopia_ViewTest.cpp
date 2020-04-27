#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_Views.hpp"
#include "utopia_DeviceReduce.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DeviceTensorContraction.hpp"

#include <utility>

namespace utopia {

    class ViewTest {
    public:
        using Scalar = double;

        ViewTest()
        {}

        void run()
        {
            UTOPIA_RUN_TEST(failing_eigen_test);
            UTOPIA_RUN_TEST(array_view_test);
            UTOPIA_RUN_TEST(static_array_view_test);
            UTOPIA_RUN_TEST(vector_view_test);
            UTOPIA_RUN_TEST(matrix_view_test);
            UTOPIA_RUN_TEST(mv_view_test);
            UTOPIA_RUN_TEST(view_assign_test);
            UTOPIA_RUN_TEST(view_binary_test);
            UTOPIA_RUN_TEST(view_unary_test);
            UTOPIA_RUN_TEST(view_composite_test);
            UTOPIA_RUN_TEST(view_inv_test);
            UTOPIA_RUN_TEST(view_transpose_test);
            UTOPIA_RUN_TEST(view_trace_test);
            UTOPIA_RUN_TEST(view_norm_test);
            UTOPIA_RUN_TEST(view_eig_2_test);
            UTOPIA_RUN_TEST(view_eig_3_test);
            UTOPIA_RUN_TEST(view_diag_test);
            UTOPIA_RUN_TEST(choose_type);
            UTOPIA_RUN_TEST(size_test);
            UTOPIA_RUN_TEST(strain_test);
            UTOPIA_RUN_TEST(inner_test);
            UTOPIA_RUN_TEST(eigen_test);
            UTOPIA_RUN_TEST(composite_expr_test);
            UTOPIA_RUN_TEST(tensor4_test);
        }

        void array_view_test()
        {
            std::vector<Scalar> a_v(2, 1.0);
            ArrayView<Scalar> a(&a_v[0], a_v.size());
            device::axpy(4.0, a, a);
        }

        void static_array_view_test()
        {
            ArrayView<Scalar, 2> a;
            ArrayView<Scalar, 2, 2> b;

            device::fill(2.0, a);
            device::fill(3.0, b);

            device::axpy(4.0, a, a);

            device::axpy(4.0, b, b);
        }

        void vector_view_test()
        {
            StaticVector2<Scalar> a;
            a.set(1.0);

            Scalar dot_a = dot(a, a);
            utopia_test_assert(approxeq(2.0, dot_a));
        }

        void matrix_view_test()
        {
            StaticMatrix<Scalar, 2, 2> a;
            a.set(1.0);

            Scalar dot_a = dot(a, a);
            utopia_test_assert(approxeq(4.0, dot_a));
        }

        void mv_view_test()
        {
            StaticMatrix<Scalar, 3, 2> A;
            StaticVector2<Scalar> x;
            StaticVector3<Scalar> y;
            StaticVector3<Scalar> expected;

            A.set(1.0);
            x.set(2.0);

            expected.set(4.0);

            y = A * x;

            utopia_test_assert(approxeq(expected, y, 1e-10));
        }

        void view_assign_test()
        {
            using V = utopia::StaticVector3<Scalar>;
            V x;
            V y;
            x.set(0.0);
            y.set(1.0);

            DeviceAssign<V, V> va(x, y);
            va.apply();

            utopia_test_assert(approxeq(x, y));
        }

        void view_binary_test()
        {
            using V = utopia::StaticVector3<Scalar>;
            V x;
            V y;
            x.set(0.0);
            y.set(1.0);

            DeviceAssign<V, DeviceBinary<V, V, Plus>> va(x, DeviceBinary<V, V, Plus>(y, y));
            va.apply();

            y.scale(2.0);
            utopia_test_assert(approxeq(x, y));
        }

        void view_unary_test()
        {
            using V = utopia::StaticVector3<Scalar>;
            V x, y;
            x.set(2.0);
            y.set(4.0);

            x = pow2(x);
            x = abs(x);

            utopia_test_assert(approxeq(x, y));
        }

        void view_composite_test()
        {
            using V3 = utopia::StaticVector3<Scalar>;
            using V2 = utopia::StaticVector2<Scalar>;
            using Mat2x3 = utopia::StaticMatrix<Scalar, 2, 3>;
            using Mat3x2 = utopia::StaticMatrix<Scalar, 3, 2>;
            using Mat2x2 = utopia::StaticMatrix<Scalar, 2, 2>;

            V3 x, y;

            x.set(2.0);
            y.set(6.0);

            x = pow2(x) + abs(x);

            utopia_test_assert(approxeq(x, y));

            Mat2x3 A;
            A.set(2.0);
            x.set(1.0);

            V2 w; w.set(1.0);
            V2 z = A * x + w;
            w.set(7.0);

            utopia_test_assert(approxeq(z, w));

            Mat3x2 A_t = transpose(A);

            Mat3x2 expected; expected.set(2.0);

            utopia_test_assert(approxeq(expected, A_t));

            A(0, 0) = 1.0;
            A(0, 1) = 2.0;
            A(0, 2) = 3.0;

            A(1, 0) = 3.0;
            A(1, 1) = 2.0;
            A(1, 2) = 1.0;

            Mat2x2 AAt = A * transpose(A);

            double det_A = det(AAt);
            utopia_test_assert(approxeq(96.0, det_A, 1e-8));
        }

        void view_transpose_test()
        {
            using Mat2x3 = utopia::StaticMatrix<Scalar, 2, 3>;
            using Mat3x2 = utopia::StaticMatrix<Scalar, 3, 2>;

            Mat2x3 A;

            A(0, 0) = 1.0;
            A(0, 1) = 2.0;
            A(0, 2) = 3.0;

            A(1, 0) = 3.0;
            A(1, 1) = 2.0;
            A(1, 2) = 1.0;


            Mat3x2 A_t = transpose(A);
            Mat3x2 A_t_expected;

            A_t_expected(0, 0) = 1.0;
            A_t_expected(1, 0) = 2.0;
            A_t_expected(2, 0) = 3.0;

            A_t_expected(0, 1) = 3.0;
            A_t_expected(1, 1) = 2.0;
            A_t_expected(2, 1) = 1.0;

            utopia_test_assert(approxeq(A_t, A_t_expected));
            utopia_test_assert(transpose(transpose(A)).is_alias(A));
        }

        void view_inv_test()
        {
            using Mat2x2 = utopia::StaticMatrix<Scalar, 2, 2>;

            Mat2x2 A;
            A(0, 0) = 1.0;
            A(0, 1) = 2.0;

            A(1, 0) = 3.0;
            A(1, 1) = 2.0;

            Mat2x2 A_inv = inv(A);
            Mat2x2 Id = A * A_inv;
            Mat2x2 expected_Id; expected_Id.set(0.0);

            expected_Id(0, 0) = 1.0;
            expected_Id(1, 1) = 1.0;

            utopia_test_assert(approxeq(Id, expected_Id, 1e-8));
            utopia_test_assert(inv(inv(A)).is_alias(A));
        }

        void view_trace_test()
        {
            using Mat2x2 = utopia::StaticMatrix<Scalar, 2, 2>;
            Mat2x2 A; A.set(2.0);
            utopia_test_assert(approxeq(trace(A), 4.0, 1e-8));
        }

        void view_norm_test()
        {
            using V2 = utopia::StaticVector2<Scalar>;
            using Mat2x3 = utopia::StaticMatrix<Scalar, 2, 3>;

            V2 v; v.set(2.0);

            double n1 = norm1(v);
            double n2 = norm2(v);
            double n_infty = norm_infty(v);

            utopia_test_assert(approxeq(4.0, n1));
            utopia_test_assert(approxeq(std::sqrt(8.0), n2, 1e-8));
            utopia_test_assert(approxeq(2.0, n_infty));

            Mat2x3 m; m.set(2.0);

            n1 = norm1(m);
            n2 = norm2(m);
            n_infty = norm_infty(m);

            utopia_test_assert(approxeq(2*3*2.0, n1));
            utopia_test_assert(approxeq(std::sqrt(2*3*4.0), n2, 1e-8));
            utopia_test_assert(approxeq(2.0, n_infty));
        }

        void view_eig_2_test()
        {
            using V2 = utopia::StaticVector2<Scalar>;
            using Mat2x2 = utopia::StaticMatrix<Scalar, 2, 2>;
            Mat2x2 A; A.set(2.0);
            A(0, 0) = 1;

            V2 e;
            Mat2x2 v;

            eig(A, e, v);

            //Oracle from MATLAB
            utopia_test_assert( approxeq(e[0], -0.561552812808830, 1e-10) );
            utopia_test_assert( approxeq(e[1],  3.561552812808830, 1e-10) );

            //first vector
            utopia_test_assert( approxeq(std::abs(v(0,0)),   0.788205438016109, 1e-10) );
            utopia_test_assert( approxeq(std::abs(v(1,0)),   0.615412209402636 , 1e-10) );

            //second vector
            utopia_test_assert( approxeq(std::abs(v(0,1)),  0.615412209402636, 1e-10) );
            utopia_test_assert( approxeq(std::abs(v(1,1)),  0.788205438016109, 1e-10) );
        }

        void view_eig_3_test()
        {
            using V3 = utopia::StaticVector3<Scalar>;
            using Mat3x3 = utopia::StaticMatrix<Scalar, 3, 3>;
            Mat3x3 A; A.set(2.0);
            A(0, 0) = 1;
            A(2, 2) = 3;

            V3 e;
            Mat3x3 v;

            eig(A, e, v);

            Mat3x3 A_actual = v * diag(e) * transpose(v);

            utopia_test_assert(approxeq(A, A_actual, 1e-8));
        }

        void view_diag_test()
        {
            using V2 = utopia::StaticVector2<Scalar>;
            using Mat2x2 = utopia::StaticMatrix<Scalar, 2, 2>;
            Mat2x2 A; A.set(2.0);

            V2 d = diag(A);

            V2 expected; expected.set(2.0);

            utopia_test_assert(approxeq(d, expected));
        }

        void choose_type()
        {
            StaticMatrix<Scalar, 2, 2> A;

            auto a_t = transpose(A);
            using AT = decltype(a_t);


            static_assert(Traits<AT>::Order == 2, "must be 2nd order tensor");

            auto expr = transpose(A) + A;
            using E = decltype(expr);

            DeviceNumber<Scalar> num;
            using T = ChooseType<DeviceNumber<Scalar>, E, E>::Type;

            static_assert(Traits<T>::Order == 2, "must be 2nd order tensor");


            auto axA = 0.5 * A;
            using AXA = decltype(axA);
            static_assert(Traits<AXA>::Order == 2, "must be 2nd order tensor");
        }

        void size_test()
        {
            StaticMatrix<Scalar, 2, 2> A;
            rows( transpose(A) );
            utopia_test_assert((rows( transpose(A) + A )) == 2);
        }

        void strain_test()
        {
            StaticMatrix<Scalar, 2, 2> A, E;
            auto expr = 0.5 * (transpose(A) + A);
            E = expr;

            static_assert(Traits<decltype(expr)>::Order == 2, "must be 2nd order tensor");
        }

        void eigen_test()
        {
            //[ 0, 0, 0;
            // 0, -0.0000088729833462074205, -0.0000044364916731037103;
            // 0, -0.0000044364916731037103, 0 ]

            StaticMatrix<Scalar, 3, 3> A, V;
            A.set(0.0);

            A(1, 1) = -0.0000088729833462074205;
            A(1, 2) = -0.0000044364916731037103;
            A(2, 1) = -0.0000044364916731037103;

            StaticVector<Scalar, 3> e;
            eig(A, e, V);

            A.set(0.0);
            A(0,0) = 0.020000000127017;
            A(0,1) = A(1,0) = 0.000000000063508;
            A(0,2) = A(2,0) = 0.000000000063508;

            eig(A, e, V);

            Scalar sum_v = sum(V);

            utopia_test_assert(sum_v == sum_v);

            A.set(0.0);
            A(0,0) = 0.02;
            A(0,1) = A(1,0) = 6.35083e-11;
            A(0,2) = A(2,0) = 6.35083e-11;

            // disp("------");
            // disp(A);

            V.set(0.0);
            eig(A, e, V);

            // disp("------");

            // disp(e);
            // disp(V);

            sum_v = sum(V);
            utopia_test_assert(sum_v == sum_v);

            utopia_test_assert( approxeq((V) * diag(e) * transpose(V), A) );
        }

        void failing_eigen_test()
        {
            StaticMatrix<Scalar, 3, 3> A, V;
            A.set(0.0);

            StaticVector<Scalar, 3> e;

            A.raw_type().copy(
            {
                0.020000000000000, 0.000001250000000, 0.000001250000000,
                0.000001250000000, 0.000004436491673, 0.0,
                0.000001250000000, 0.0, 0.000004436491673
            });

            eig(A, e, V);
        }

        void inner_test()
        {
            StaticMatrix<Scalar, 2, 2> A, E;
            auto expr = 0.5 * (transpose(A) + A);
            utopia_test_assert( rows(expr) == 2 );
        }

        void composite_expr_test()
        {
            StaticMatrix<Scalar, 2, 2> F, S;
            F.set(0.0);
            S.set(0.0);

            F(0,0) = -0.5;
            F(0,1) = -0.5;

            S(0,0) = 5.5999999999999996;
            S(1,1) = 2.3999999999999999;

            const Scalar SdotE = inner(S, 0.5 * (F + transpose(F)));
            utopia_test_assert(approxeq(-2.8000, SdotE, 1e-8));
        }

        void tensor4_test()
        {
            StaticMatrix3x3<Scalar> m1, m2, m3;
            m1.identity();
            m2.identity();
            m3.identity();

            //or Tensor4th<Scalar, 3, 3, 3, 3>
            Tensor3x3x3x3<Scalar> t;
            t.identity();

            //get/set
            t(0, 1, 2, 0) = 120;

            t = t + 0.5 * t;
            t = t - 1.5 * t;

            //t_{ijkl} = m_{ij} * m_{kl}
            t = tensor_product<0, 1, 2, 3>(m1, m1);
            // disp(t);

            // t_{ijkl} = m_{ik} * m_{jl} ==
            // t_{ikjl} = m_{ij} * m_{kl}
            t = tensor_product<0, 2, 1, 3>(m1, m1);

            // m2_{ij} = t_{ijkl} * m_{kl}
            m2 = contraction(t, m1);
            // disp(m2);

            Scalar val = inner(m1, contraction(tensor_product<0, 2, 1, 3>(m1, m1), m1));
            // disp(val);

            Tensor3x3x3x3<Scalar> t_mult = t * t;


            static const int Dim = 2;
            Tensor4th<Scalar, Dim, Dim, Dim, Dim> C, Id, actual;
            C.set(0.0);
            Id.identity();

            auto kronecker_delta = [](const SizeType &i, const SizeType &j) -> bool
            {
                return (i==j) ? 1.0 : 0.0;
            };

            for(SizeType i = 0; i < Dim; ++i) {
                for(SizeType j = 0; j < Dim; ++j) {
                    for(SizeType k = 0; k < Dim; ++k) {
                        for(SizeType l = 0; l < Dim; ++l) {
                            Scalar val = 120 * kronecker_delta(i,j)* kronecker_delta(k,l);
                            val += 80 * (kronecker_delta(i,k)* kronecker_delta(j,l));
                            val += 80 * (kronecker_delta(i,l)* kronecker_delta(j,k));
                            C(i, j, k, l) = val;
                        }
                    }
                }
            }

            actual = Id * C;

            // disp(actual);
            // disp("----------");
            // disp(C);

            utopia_test_assert(approxeq(actual, C, 1e-10));
        }


    };

    void view()
    {
        ViewTest().run();
    }

    UTOPIA_REGISTER_TEST_FUNCTION(view);
}

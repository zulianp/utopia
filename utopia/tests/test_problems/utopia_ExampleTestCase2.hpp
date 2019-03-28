/*! \file utopia_ExampleTestCase2.hpp
     Nonlinear Semismooth Newton method
     Created by Alena Kopanicakova
*/

#ifndef UTOPIA_SOLVER_EXAMPLE_TEST_CASE2_HPP
#define UTOPIA_SOLVER_EXAMPLE_TEST_CASE2_HPP

#include "utopia_Core.hpp"
#include <vector>


namespace utopia
{
    /**
     * @brief      Example used for testing Nonlinear Semismooth method.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class ExampleTestCase2
    {
        typedef typename utopia::Traits<Vector>::Scalar Scalar;

    public:

        ExampleTestCase2()
        {
        }

        /**
         * @brief      Assembles operators.
         *
         * @param[in]  N     The problem size.
         */
        void assembly(SizeType N)
        {
            const Scalar h = 1.0 / (N - 1);
            A = sparse(N, N, 3);

            {
                Write<Matrix> w(A);
                Range rr = row_range(A);

                for (SizeType i = rr.begin(); i != rr.end(); i++) {
                    if(i == 0 || i == N - 1) continue;

                    const Scalar inv_2h = (1. / (h * h));

                    // diag
                    A.set(i, i, 2.0 * inv_2h);

                    // upper diag
                    A.set(i, i + 1, -1.0 * inv_2h);

                    // lower diag
                    A.set(i, i - 1, -1.0 * inv_2h);
                }

                if(rr.begin() == 0) {
                    A.set(0, 0, 1.0);
                }

                if(rr.inside(N - 1)) {
                    A.set(N - 1, N - 1, 1.0);
                }
            }

            g = zeros(N);
            Vector s = zeros(N);
            Range rr = range(s);
            {
                Write<Vector> w(s);
                Scalar a = 0.;
                for (SizeType i = rr.begin() + 1; i != rr.end(); i++) {
                    a += h;
                    s.set(i, a);
                }
            }

            {
                Write<Vector> w(g);
                Read<Vector> r(s);

                Range rg = range(g);
                for (SizeType i = rg.begin(); i != rg.end(); i++) {
                    Scalar ss = 0.5 + ((s.get(i) - 0.5) * (s.get(i) - 0.5));
                    g.set(i, ss);
                }
            }

            b = values(N, 10.0);
            {
                Write<Vector> w(b);

                Range b_range = range(b);
                if(b_range.begin() == 0) b.set(0, 0);
                if(b_range.end() == N) b.set(N - 1, 0);
            }
        }

        /**
         * @brief      Gets the operators.
         *
         * @param[in]  N      The problem size.
         * @param      A_out  The stifness matrix.
         * @param      b_out  The constraints.
         * @param      g_out  The upper bound.
         */
        void getOperators(const SizeType N, Matrix &A_out, Vector &b_out, Vector &g_out) {
            this->assembly(N);
            A_out = A;
            b_out = b;
            g_out = g;
        }


    private:
        Matrix A; /*!< Laplacian  */
        Vector b;
        Vector g; /*!< Gradient  */
    };

}
#endif // UTOPIA_SOLVER_EXAMPLE_TEST_CASE2_HPP

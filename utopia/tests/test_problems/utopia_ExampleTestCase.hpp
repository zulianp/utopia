/*! \file utopia_ExampleTestCase.hpp
     Nonlinear Semismooth Newton method 
     Created by Alena Kopanicakova 
*/

#ifndef UTOPIA_SOLVER_EXAMPLE_TEST_CASE_HPP
#define UTOPIA_SOLVER_EXAMPLE_TEST_CASE_HPP

#include "utopia_Core.hpp"
#include <vector>


namespace utopia 
{
    /**
     * @brief      Example of assembly 2D Laplace problems. 
     *
     * @tparam     Matrix 
     * @tparam     Vector 
     */
    template<class Matrix, class Vector>
    class ExampleTestCase  
    {
        typedef typename utopia::Traits<Vector>::Scalar Scalar;
        typedef typename utopia::Traits<Vector>::SizeType SizeType;

    public:
        ExampleTestCase()
        {
        }
      
      /**
       * @brief      Assembles 
       *
       * @param[in]  N     The problem size. 
       */
      void assembly(SizeType N)
      {

            Scalar h = 1.0 / (N - 1);
            A = sparse(N, N, 3);
            {
                Write<Matrix> w(A);
                Range rr = row_range(A);

                
                for (SizeType i = rr.begin(); i != rr.end(); i++) {
                    const SizeType ip1 = i+1;
                    const Scalar inv2h = (1 / (h * h));

                    // diag 
                    A.set(i, i, 2.0 * inv2h);

                    // upper diag
                    if(ip1 < N) {
                        A.set(i, i + 1, -1.0 * inv2h);
                    }

                    // lower diag
                    if (ip1 < rr.end()) {
                        A.set(i + 1, i, -1.0 * inv2h);
                    }
                }

                if(rr.begin() == 0) {
                    A.set(0, 0, 1.0);
                    for (SizeType i = 0; i != N; i++) {
                        A.set(0, i, 0);
                    }
                }

                if(rr.inside(N-1)) {
                    A.set(N - 1, N - 1, 1.0);
                    for (SizeType i = 0; i != N; i++) {
                        A.set(N - 1, i, 0);
                    }
                }
            }

            B = identity(N, N);
            {
                Write<Matrix> w (B);
                Range B_range = row_range(B);
                if(B_range.begin() == 0)    B.set(0,0, 0); 
                if(B_range.end() == N)    B.set(N - 1, N - 1, 0);
            }


            Scalar periods = 4;
            upbo = zeros(N);
            {
                Write<Vector> w (upbo);
                Scalar help_x = 0.0;
                Range upbo_range = range(upbo);
                for (SizeType i = upbo_range.begin(); i != upbo_range.end() ; i++){
                    help_x += h;
                    Scalar c = 2.0 * 3.14 * periods * (2.0*help_x - 1.0)/2.0;
                    Scalar s = 0.5 + ((std::cos(c)) - 0.5) * ((std::cos(c)) - 0.5);
                    upbo.set(i, s);
                }
            }
        }

    

    /**
     * @brief      Gets the operators ( )
     *
     * @param[in]  N         The problem size. 
     * @param      A_out     The stifness matrix.
     * @param      B_out     The 
     * @param      upbo_out  The upper bound.
     */
    void getOperators(const SizeType N, Matrix &A_out, Matrix &B_out, Vector &upbo_out)
    {
        this->assembly(N);
        A_out = A;
        B_out = B;
        upbo_out = upbo;
    }


    private:
        Matrix A; /*!< Laplacian  */  
        Matrix B; /*!< Boundary operator  */  
        Vector upbo; /*!< Constraints  */  
    };

}
#endif // UTOPIA_SOLVER_EXAMPLE_TEST_CASE_HPP

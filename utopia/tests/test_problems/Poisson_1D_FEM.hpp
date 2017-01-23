/*
* @Author: alenakopanicakova
* @Date:   2016-03-25
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-07-29
*/
#ifndef UTOPIA_1D_FEM_EXAMPLE_HPP
#define UTOPIA_1D_FEM_EXAMPLE_HPP

#include <vector>
#include "utopia_Function.hpp"




namespace utopia 
{
    /**
     * @brief      Example of assembly of 1D FEM Laplacian. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class Poisson_1D  
    {
        typedef typename utopia::Traits<Vector>::Scalar Scalar;
        typedef typename utopia::Traits<Vector>::SizeType SizeType;

    public:
        Poisson_1D()
        {
            
        }

        /**
         * @brief      Assembles stifness matrix and right hand side. 
         *
         * @param[in]  N     The size of problem. 
         *
         * @return     
         */
        bool assemble(const SizeType & N)
        {
            get_L(N); 
            get_rhs(N);
            return true;   
        }

         /**
          * @brief      Gets the stifness matrix
          *
          * @param[in]  N     THe size of problem.
          */
        void get_L(const SizeType & N)
        {
            Scalar h = 1.0 / (N - 1);
            A = sparse(N, N, 3);
            {
                Write<Matrix> w(A);
                Range rr = rowRange(A);
                Range cr = colRange(A);
                for (SizeType i = rr.begin(); i != rr.end(); i++) 
                {
                    const SizeType ip1 = i+1;
                    const Scalar inv2h = (1 / (h * h));

                    // diag 
                    A.set(i, i, 2.0 * inv2h);

                    // upper diag
                    if(ip1 < cr.end()) {
                        A.set(i, i + 1, -1.0 * inv2h);
                    }

                    // lower diag
                    if (ip1 < rr.end()) {
                        A.set(i + 1, i, -1.0 * inv2h);
                    }
                }
            }
        }


        /**
         * @brief      Get the Mass matrix. 
         *
         * @param[in]  N     THe problem size.
         * @param      M     The mass matrix. 
         */
        void get_M(const SizeType & N, Matrix &M)
        {
            Scalar h = 1.0 / (N - 1);
            M = sparse(N, N, 3);
            {
                Write<Matrix> w(M);
                Range rr = rowRange(M);
                Range cr = colRange(M);
                for (SizeType i = rr.begin(); i != rr.end(); i++) 
                {
                    const SizeType ip1 = i + 1;
                    const Scalar h6 = (h / 6);

                    // diag 
                    M.set(i, i, 4.0 * h6);

                    // upper diag
                    if(ip1 < cr.end()) {
                        M.set(i, i + 1, 1.0 * h6);
                    }

                    // lower diag
                    if (ip1 < rr.end()) {
                        M.set(i + 1, i, 1.0 * h6);
                    }
                }
            }
        }
    



        /**
         * @brief      Gets the rhs.
         *
         * @param[in]  N     The problem size.
         */
        void get_rhs(const SizeType & N)
        {
            Scalar h = 1.0 / (N - 1);
            rhs = zeros(N);
            {
                Write<Vector> w (rhs);
                Range rhs_range = range(rhs);
                Scalar x_step = 0.0, source;
                for (SizeType i = rhs_range.begin(); i != rhs_range.end() ; i++)
                {
                    x_step += h;
                    source = std::sin(3.16 * x_step);
                    rhs.set(i, source);
                }
            }

            get_M(N, M); 
            rhs = M * rhs; 
        }



        /**
         * @brief      Get the operators.
         *
         * @param[in]  N     { problem size }
         * @param      L     { stiffness }
         * @param      rhs   { right hand side  }
         */
        void getOperators(const SizeType & N_in, Matrix & A_in, Vector &rhs_in)
        {
            this->assemble(N_in);
            A_in    =   A;
            rhs_in  =   rhs; 
        }


        /**
         * @brief      Gets the transfer operators.
         *
         * @param[in]  N_in  The size of problem.
         * @param      R_in  The restriction operator. 
         * @param      I_in  The interpolation operator. 
         *
         * @return     The transfer operators.
         */
        bool getTransferOperators(const SizeType & N_in, Matrix & R_in, Matrix & I_in)
        {
            assemble_transfer(N_in); 
            R_in = R; 
            I_in = I; 
            return true; 
        }


        /**
         * @brief      Assembles  the interpolation/restriction for 1D grid.
         *
         * @param[in]  N    The size of problems. 
         *
         * @return  
         */
        bool assemble_transfer(const SizeType & N)
        {

            I = sparse(N, (N - 1)/2, 4);
            {
                Write<Matrix> w(I);
          //      Range rr = rowRange(I);
                Range cr = colRange(I);
                SizeType ip1 = 0; 

                for (SizeType i = cr.begin(); i != cr.end(); i++) 
                {
                    I.set(ip1, i, 1);
                    I.set(ip1 + 1, i, 2);
                    I.set(ip1 + 2, i, 1);
                    ip1 = i + 3;
                }
            }

            disp(I); 

            I = 0.5 * I; 
            R = transpose(I);
            R *= 0.25; 

            return true; 
        }


    private:
        Matrix A;       /*!< stiffness  matrix     */  
        Matrix M;       /*!< mass matrix           */  
        Vector rhs;     /*!< rhs                   */  


        Matrix R;       /* restriction              */ 
        Matrix I;       /* interpolation            */
        
    };

}
#endif // UTOPIA_1D_FEM_EXAMPLE_HPP

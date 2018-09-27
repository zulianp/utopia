#ifndef UTOPIA_LBFGS_HPP
#define UTOPIA_LBFGS_HPP

#include "utopia_Core.hpp"


namespace utopia
{

    template<class SparseMatrix, class DenseMatrix, class Vector>
    class LBFGS : public HessianApproximation<SparseMatrix, Vector>
    {

        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        public:

            LBFGS(const SizeType & m): m_(m), current_m_(0)
            {

            }

            virtual bool initialize(Function<SparseMatrix, Vector> &fun, const Vector &x, SparseMatrix &H) override
            {
                SizeType n = local_size(x).get(0);

                //for sparse
                if(!fun.hessian(x, H))
                    H = local_identity(n, n);

                current_m_ = 0; 

                return true;
            }


                virtual bool approximate_hessian(
                    Function<SparseMatrix, Vector> &fun,
                    const Vector &sol_new,
                    const Vector &step,
                    SparseMatrix &hessian_old_new,
                    Vector &grad_old_new) override
                {


                    return false; 
                }




            void set_memory_size(const SizeType & m)
            {
                m_ = m; 
            }

            SizeType get_memory_size()
            {
                return m_; 
            }


        // void petsc_matrix_permuation()
        // {

        // auto k = 15;
        // auto m = 4;

        // DMatrixd A = local_values(k, m, 2.);
        // DMatrixd A_per=local_values(k, m, 0.); 


        // {
        //     Write<DMatrixd>  w(A); 

        //     A.set(0,0., 1); 
        //     A.set(1,1. ,3); 
        //     A.set(2,2. ,7); 
        //     A.set(11,3., 5); 
        //     A.set(13,0., 9); 
        // }


        // // shift rows...
        // {
        //     Write<DMatrixd>  write(A_per); 
        //     Read<DMatrixd>  read(A); 

        //     Range r = row_range(A);
        //     Range c = col_range(A);

        //     for(SizeType i = r.begin(); i != r.end(); ++i)
        //     {
        //         for(SizeType j = c.begin(); j != c.end(); ++j)
        //         {
        //             if(j < m-1)
        //                 A_per.set(i, j, A.get(i, j+1)); 
        //         }
        //     }
        // }

        // DVectord v = local_values(k, 999); 

        // // set last row to be new vector 
        // {
        //     Write<DMatrixd>  write(A_per); 
        //     Read<DVectord>  read(v); 

        //     Range r = row_range(A);
        //     Range c = col_range(A);

        //     for(SizeType i = r.begin(); i != r.end(); ++i)
        //     {
        //         A_per.set(i, c.end()-1, v.get(i)); 
        //     }
        // }

        // disp(A_per); 


        // DMatrixd AtA = transpose(A_per)*A_per; 

        // disp(AtA); 
    
        // }


        private:
            static_assert(utopia::is_sparse<SparseMatrix>::value, "BFGS does not support sparse matrices."); 
            static_assert(!utopia::is_sparse<DenseMatrix>::value, "BFGS does not support sparse matrices."); 

            SizeType m_; 
            SizeType current_m_; 

        };



}

#endif //UTOPIA_LBFGS_HPP

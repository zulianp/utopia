#ifndef UTOPIA_LBFGSB_HPP
#define UTOPIA_LBFGSB_HPP

#include "utopia_Core.hpp"


namespace utopia
{

template<class Matrix, class Vector>
class LBFGSB : public HessianApproximation<Matrix, Vector>
{

    typedef UTOPIA_SCALAR(Vector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

    public:

        LBFGSB(const SizeType & m): m_(m), current_m_(0)
        {

        }

        virtual bool initialize(Function<Matrix, Vector> &fun, const Vector &x) override
        {
            SizeType n = local_size(x).get(0);
            H0_ = local_identity(n, n);

            // initialize some other vectors and matrices

            this->initialized(true);
            current_m_ = 0; 

            return true;
        }

        virtual bool update(const Vector &  s, const Vector &  y) override
        {

            if(!this->initialized())
            {
                utopia_error("BFGS::update: Initialization needs to be done before updating. \n"); 
                return false; 
            }

            theta_ = dot(y,y)/dot(y,s); 


            // this needs to be done durign first it, in order to initialize matrices 
            if(current_m_ == 0)
            {
                this->init_mat_from_vec(Y_, y); 
                this->init_mat_from_vec(S_, s); 
            }
            // here we are just building first m elements of matrices
            else if(current_m_ < m_)
            {
                this->add_col_to_mat(Y_, y); 
                this->add_col_to_mat(S_, s); 
            }
            // here, we are removing first col, shifting, appending new 
            else
            {  
                this->shift_cols_left_replace_last_col(Y_, y); 
                this->shift_cols_left_replace_last_col(S_, s); 

            }




            current_m_++; 

            return true;
        }

        virtual bool apply_Hinv(const Vector & /* g */, Vector & /*s */) override
        {
            // TODO
            return true;
        }


        virtual Matrix & get_Hessian() override
        {
            std::cerr << "--- not implemented yet---- \n";
            return H0_;
        }

        virtual Matrix & get_Hessian_inv() override
        {
            std::cerr << "--- not implemented yet---- \n";
            return H0_;
        }


        void set_memory_size(const SizeType & m)
        {
            m_ = m;
        }

        SizeType get_memory_size()
        {
            return m_;
        }


    private:
        /**
         * @brief Appends vector into matrix as new last column
         * 
         * @param M Matrix to be appended to
         * @param col Vector to be added
         */
        void add_col_to_mat(Matrix & M, const Vector & col)
        {
            const Matrix M_old = M;
            M = values(size(M_old).get(0), size(M_old).get(1) + 1, 0.0);

            {   // begin lock

                Write<Matrix>  write(M);

                Read<Matrix>   read(M_old);
                Read<Vector>      read_vec(col);

                Range r = row_range(M_old);
                Range c = col_range(M_old);

                Range r_new = row_range(M);
                Range c_new = col_range(M);

                for (SizeType i = r.begin(); i != r.end(); ++i) {
                    for (SizeType j = c.begin(); j != c.end(); ++j)
                        M.set(i, j, M_old.get(i, j));
                }

                for (SizeType i = r_new.begin(); i != r_new.end(); ++i)
                    M.set(i, c_new.end() - 1, col.get(i));

            } // end of lock
        }

        /**
         * @brief Add shift columns of the matrix to the left by one col. 
         * Then, replace last column with new vector
         * 
         * @param M Matrix to be shifted
         * @param col Vector to be added to the last col
         */
        void shift_cols_left_replace_last_col(Matrix & M, const Vector & col)
        {
            const Matrix M_old = M;

            {   // begin lock

                Write<Matrix>  write(M);

                Read<Matrix>   read(M_old);
                Read<Vector>      read_vec(col);

                Range r = row_range(M_old);
                Range c = col_range(M_old);

                for (SizeType i = r.begin(); i != r.end(); ++i) {
                    for (SizeType j = c.begin(); j != c.end() - 1; ++j) {
                        M.set(i, j, M_old.get(i, j + 1));
                    }
                    M.set(i, c.end() - 1, col.get(i));
                }

            } // end of lock
        }


        void init_mat_from_vec(Matrix & M, const Vector & col)
        {
            M = values(size(col).get(0), 1, 0.0);

            {   // begin lock

                Write<Matrix>  write(M);
                Read<Vector>   read_vec(col);

                Range r = row_range(M);
                Range c = col_range(M);

                for (SizeType i = r.begin(); i != r.end(); ++i)
                    M.set(i, c.end()-1, col.get(i));

            } // end of lock
        }


    private:
        static_assert(utopia::is_sparse<Matrix>::value, "BFGS does not support sparse matrices.");

        SizeType m_; // memory size
        SizeType current_m_; // current amount of vectors in the memory
        Matrix H0_;  // identity

        Scalar theta_; // scaling param for H0

        Matrix Y_;  // matrix n \times m of grad diffs, Y_k = [y_{k-m}, ..., y_{k-1}], where y_k = g_{k+1} - g_k
        Matrix S_;  // matrix n \times m of directions, S_k = [s_{k-m}, ..., s_{k-1}], where s_k = x_{k+1} - x_k
        
        Matrix W_;  // W_k = [Y_k \theta S_k]
        Matrix M_;  /* M_k  =   | -D       L_k^T  |
                                | L_k      \theta S_k^T S_k | */

};



}

#endif //UTOPIA_LBFGSB_HPP

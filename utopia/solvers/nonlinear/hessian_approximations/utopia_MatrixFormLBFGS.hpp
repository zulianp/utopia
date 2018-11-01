#ifndef UTOPIA_MatrixFormLBFGS_HPP
#define UTOPIA_MatrixFormLBFGS_HPP

#include "utopia_Core.hpp"

namespace utopia
{

template<class Matrix, class DenseMatrix, class Vector>
class MatrixFormLBFGS : public HessianApproximation<Vector>
{

    typedef UTOPIA_SCALAR(Vector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

    typedef utopia::LinearSolver<Matrix, Vector>    LinSolver;

    public:

        MatrixFormLBFGS( const SizeType & m, 
                const std::shared_ptr <LinSolver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >()):
                m_(m), current_m_(0), linear_solver_(linear_solver)
        {
           
        }


        virtual void initialize(const SizeType & n) override
        {
            H_ = local_identity(n, n);

            current_m_ = 0; 

            SizeType local_col_W = 0; 
            if(mpi_world_rank()==0)
                local_col_W = 2.0 * m_; 

            // W \in R^{n \times 2m }
            W_ = local_values(n, local_col_W, 0.0); 

            // M \in R^{2m \times 2m }
            M_ = local_values(local_col_W, local_col_W, 0.0); 

            SizeType local_col_YS = 0; 
            if(mpi_world_rank()==0)
                local_col_YS = m_; 

            Y_ = local_values(n, local_col_YS, 0.0); 
            S_ = local_values(n, local_col_YS, 0.0); 

            theta_ = 1.0; 

            d_elements_.resize(m_); 
            L_dots_.resize(m_, std::vector<Scalar>(m_));

            this->initialized(true);
        }


        virtual void set_linear_solver(const std::shared_ptr<LinSolver> &linear_solver)
        {
            linear_solver_ = linear_solver; 
        }


        virtual bool update(const Vector &  s, const Vector &  y) override
        {
            if(!this->initialized())
            {
                utopia_error("BFGS::update: Initialization needs to be done before updating. \n"); 
                return false; 
            }

            Scalar denom    = dot(y,s); 
            Scalar nom      = dot(y,y);

            // if denom > eps, hessian approx. should be positive semidefinite
            if(denom < 1e-10)
            {
                utopia_warning("L-BFGS-B: Curvature condition not satified. Skipping update. \n"); 
                return false; 
            }

            theta_ = nom/denom; 

            this->shift_cols_left_replace_last_col(Y_, y); 
            this->shift_cols_left_replace_last_col(S_, s); 

            this->buildL_d_elems(s); 

            this->buildW();     
            this->buildM(); 

            current_m_++; 

            return true;
        }

        virtual bool apply_Hinv(const Vector & g, Vector & s) const override
        {
            this->apply_inverse_to_vec(g, W_, s); 
            return true;
        }


        virtual bool apply_reduced_Hinv(const Vector & feasible_set, const Vector & g, Vector & s) const
        {
            // SizeType local_feasible_size = reduced_primal_method_.get_local_size_feasible_set(feasible_set); 
            // Vector g_reduced = local_zeros(local_feasible_size);         
            // reduced_primal_method_.build_reduced_vector(g, feasible_set, g_reduced); 

            // DenseMatrix W_reduced = local_values(local_feasible_size, local_size(W_).get(1), 0.0);
            // reduced_primal_method_.build_reduced_matrix(W_, feasible_set, W_reduced); 

            // this->apply_inverse_to_vec(g_reduced, W_reduced, s); 

            return true;
        }


        virtual bool apply_H(const Vector & v, Vector & result) const override
        {; 
            if(current_m_ > m_)
            {
                Vector p = transpose(W_) *v; 
                this->apply_M(p, result); 

                result = (theta_ * v) - W_*result; 
            }
            else
                result = v; 

            return false; 
        }


        virtual Matrix & get_Hessian() 
        {
            if(current_m_ > m_)
            {
                H_ = local_identity(local_size(H_).get(0), local_size(H_).get(1)); 
                H_ = theta_ * H_; 

                DenseMatrix result; 
                DenseMatrix WT = transpose(W_); 

                this->apply_M(WT, result);  
                H_ =  H_ - (W_ * result);   
            }
            else
            {
                H_ = local_identity(local_size(H_).get(0), local_size(H_).get(1)); 
            }

            return H_; 
        }


        void set_memory_size(const SizeType & m)
        {
            m_ = m;
        }

        SizeType get_memory_size() const 
        {
            return m_;
        }


    private:
        void shift_cols_left_replace_last_col(DenseMatrix & M, const Vector & col) const
        {
            const DenseMatrix M_old = M;

            {  
                Write<DenseMatrix>  write(M);

                Read<DenseMatrix>   read(M_old);
                Read<Vector>      read_vec(col);

                Range r = row_range(M_old);
                Range c = col_range(M_old);

                for (SizeType i = r.begin(); i != r.end(); ++i) {
                    for (SizeType j = c.begin(); j != c.end() - 1; ++j) {
                        M.set(i, j, M_old.get(i, j + 1));
                    }
                    M.set(i, c.end() - 1, col.get(i));
                }
            } 
        }

        void buildW()
        {
            DenseMatrix S_theta = theta_ * S_; 
            W_ = DenseMatrix(Blocks<DenseMatrix>(1, 2,
            {
                make_ref(Y_), make_ref(S_theta)
            }));
        }


        void buildL_d_elems(const Vector &s)
        {
            Vector dots = transpose(Y_) * s; 
            
            // moving things one row up and deleting diagonal elements
            for(auto i =0; i < m_-1; i++)
            {
                for(auto j =0; j < i; j++)
                {
                    L_dots_[i][j] = L_dots_[i+1][j]; 
                }
                d_elements_[i] = d_elements_[i+1]; 
            }

            {
                Read<Vector> rw(dots); 
                auto r = range(dots); 
                
                for(auto i=r.begin(); i != r.end(); ++i)
                {
                    if(i==r.end()-1)
                        d_elements_[m_-1] = dots.get(i); 
                    else
                        L_dots_[m_-1][i] = dots.get(i); 
                }
            }
        }


        void buildM()
        {           
            // transpose(S_) * S_ doesn't work with dense matrices 
            DenseMatrix SS = theta_ * transpose(S_) * S_; 

            {
                Write<DenseMatrix> wM(M_); 
                Read<DenseMatrix> rS(SS); 

                auto r = row_range(M_); 
                auto c = col_range(M_); 

                for(auto i=r.begin(); i!=r.end(); ++i)
                {
                    for(auto j=c.begin(); j!=c.end(); ++j)
                    {
                        Scalar value = 0.0; 
                        
                        if(i < m_ && j < m_ && i==j) // top left block
                        {
                            value =  -1.0 * d_elements_[j]; 
                        }
                        else if(i >= m_ && j < m_ ) // bottom left  block 
                        {
                            value = L_dots_[i -m_][j]; 
                        }
                        else if(i < m_ && j >= m_ ) // top right  block 
                        {
                            value = L_dots_[j-m_][i]; 
                        }
                        else if(i >= m_ && j >= m_ ) // bottom right  block 
                        {
                            value = SS.get(i - m_, j - m_); 
                        }

                        M_.set(i,j, value); 
                    }
                }
            }

        }


        // Sherman - Morison - Woodbury formula for computing inverse
        // this formula could be spped-up by using M_inv instead of M 
        virtual void apply_inverse_to_vec(const Vector &  g, const DenseMatrix & W, Vector & s) const 
        {
            if(current_m_ > m_)
            {
                Vector WTg = transpose(W) * g; 
                Vector MWTg; 
                this->apply_M(WTg, MWTg); 
      
                DenseMatrix WTW  = 1.0/theta_ * transpose(W) * W; 

                DenseMatrix MWTW; 
                this->apply_M(WTW, MWTW); 

                DenseMatrix N = DenseMatrix(local_identity(local_size(MWTW))) - MWTW; 
                Vector N_inv_v = local_zeros(local_size(WTg)); 

                linear_solver_->solve(N, MWTg, N_inv_v);                  

                s = (1.0/theta_ )* g; 
                s += 1.0/(theta_*theta_) * (W* N_inv_v);    
            }
            else
            {
                s = g; 
            }
        }


        void apply_M(const Vector & v, Vector & result) const 
        {
            linear_solver_->solve(M_, v, result);
        }


        void apply_M(const DenseMatrix & RHS, DenseMatrix & result) const 
        {
            result = local_values(local_size(RHS).get(0), local_size(RHS).get(1), 0.0); 
            MatLinearSolver<DenseMatrix, DenseMatrix, Vector> mat_solver(linear_solver_); 
            mat_solver.solve(M_, RHS, result); 
        }


    private:
        // static_assert(!utopia::is_sparse<Matrix>::value, "LBFGS does not support sparse matrices.");

        SizeType m_; // memory size

        std::vector<Scalar> d_elements_; 
        std::vector<std::vector<Scalar> > L_dots_; 


        SizeType current_m_; // iteration number 
        Matrix H_;  // identity

        Scalar theta_; // scaling param for H0

        DenseMatrix Y_;  // matrix n \times m of grad diffs, Y_k = [y_{k-m}, ..., y_{k-1}], where y_k = g_{k+1} - g_k
        DenseMatrix S_;  // matrix n \times m of directions, S_k = [s_{k-m}, ..., s_{k-1}], where s_k = x_{k+1} - x_k
        
        DenseMatrix W_;  // W_k = [Y_k \theta S_k]
        DenseMatrix M_;  /* M_k  =   | -D       L_k^T  |
                                | L_k      \theta S_k^T S_k | */   // in contrast to the paper, this is not the inverse

        std::shared_ptr<LinSolver> linear_solver_;     /*!< Linear solver parameters. - Ideally LU  */  

};



}

#endif //UTOPIA_MatrixFormLBFGS_HPP

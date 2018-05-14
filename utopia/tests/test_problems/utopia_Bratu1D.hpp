#ifndef UTOPIA_BRATU_1D_FD_HPP
#define UTOPIA_BRATU_1D_FD_HPP

#include "utopia_Core.hpp"
#include <vector>


namespace utopia 
{
    
    template<class Matrix, class Vector>
    class Bratu1D  : public Function<Matrix, Vector>
    {
        typedef typename utopia::Traits<Vector>::Scalar Scalar;
        typedef typename utopia::Traits<Vector>::SizeType SizeType;

    public:
        Bratu1D(const Scalar & n, const Scalar & lambda = 10): 
          n_(n), lambda_(lambda)
        {
          h_ = 1.0/(n_ - 1.0); 

          assemble_laplacian_1D(); 
        }

        bool value(const Vector &x, typename Vector::Scalar &energy) const override 
        {
          energy = 0.5 * dot(x, A_*x); 
          energy +=  lambda_ * sum(exp(x)); 

          return true;
        }

        bool gradient(const Vector &x, Vector &gradient) const override 
        {
          gradient = (A_ * x) + (lambda_ * exp(x)); 

          // enforce BC conditions
          {
            Write<Vector> w(gradient);
            Range r = range(gradient);

            if(r.begin() == 0) 
              gradient.set(0, 0.0);

            if(r.end() == n_) 
              gradient.set(n_-1, 0.0);            
          }

          return true;
        }



        bool hessian(const Vector &x, Matrix &hessian) const override 
        {
          hessian = diag(lambda_ * exp(x)); 
          hessian += A_;

          // enforce BC conditions
          {
            Range r = row_range(hessian);
            Write<Matrix> w(hessian);

            if(r.begin() == 0) 
            {
              hessian.set(0, 0, 1.0);
              hessian.set(0, 1, 0.0);
            }

            if(r.end() == n_) 
            {
              hessian.set(n_-1, n_-1, 1.0);
              hessian.set(n_-1, n_-2, 0.0);
            }
          }

          return true;        
        }


      void apply_bc_to_initial_guess(Vector &x)
      {
          // here, we assume that BC are 0 on the initial and end part of subdomain
          {
            Write<Vector> w(x);
            Range r = range(x);

            if(r.begin() == 0) 
              x.set(0, 0.0);

            if(r.end() == n_) 
              x.set(n_-1, 0.0);            
          }
      }

    private: 

        void assemble_laplacian_1D()
        {

          A_ = sparse(n_, n_, 3); 

          {
            Write<Matrix> w(A_);
            Range r = row_range(A_);

            for(SizeType i = r.begin(); i != r.end(); ++i) 
            {
              if(i > 0) 
                A_.add(i, i - 1, -1.0);    

              if(i < n_-1) 
                A_.add(i, i + 1, -1.0);


              A_.add(i, i, 2.0);
            }
          }

          A_ *= 1.0/(h_*h_); 
        }


  
    private:
        Scalar n_;  // global size
        Scalar h_; // grid size 
        Scalar lambda_; // combustion factor
        Matrix A_; // constant part of eq... 

    };

}
#endif // UTOPIA_BRATU_1D_FD_HPP

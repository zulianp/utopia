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
        Bratu1D(const Scalar & n, const Scalar & lambda = 5): 
          n_(n), lambda_(lambda)
        {
        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override 
        {

            return true;
        }

        bool gradient(const Vector &point, Vector &gradient) const override 
        {

            return true;
        }

        bool hessian(const Vector &point, Matrix &hessian) const override 
        {

            return true;        
        }

  
    private:
        Scalar n_; 
        Scalar lambda_; 
        Matrix A_; 

    };

}
#endif // UTOPIA_BRATU_1D_FD_HPP

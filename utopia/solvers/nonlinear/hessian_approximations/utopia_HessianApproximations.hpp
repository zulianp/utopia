#ifndef UTOPIA_HESSIAN_APPROXIMATIONS_HPP
#define UTOPIA_HESSIAN_APPROXIMATIONS_HPP

#include "utopia_Core.hpp"


namespace utopia
{

template<class Matrix, class Vector>
class HessianApproximation
{
    typedef UTOPIA_SCALAR(Vector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


public:

    HessianApproximation(): num_tol_(1e-12), initialized_(false)
    {

    }

    virtual ~HessianApproximation() { }

    virtual bool initialize(Function<Matrix, Vector> &fun, const Vector &x) = 0;

    Scalar num_tol()const 
    {
        return num_tol_; 
    }

    void num_tol(Scalar & tol ) 
    {
        num_tol_ = tol; 
    }

    void initialized(const bool init) 
    {
        initialized_ = init; 
    }

    bool initialized() const 
    {
        return initialized_; 
    }

    // applications of inverse of Hessian
    virtual bool apply_Hinv(const Vector & /* g */, Vector & /*s */) const  = 0;
    virtual Scalar compute_uHinvv_dot(const Vector &/*u*/, const Vector & /*v*/) const {return false; }

    // applications of Hessian
    virtual bool apply_H(const Vector & /*v*/ , Vector & /*r */) const  {return false; }
    virtual Scalar compute_uHv_dot(const Vector &/*u*/, const Vector & /*v*/) const {return false; }
    virtual Scalar compute_uHu_dot(const Vector &/*u*/) const {return false; }

    // refresh vectors
    virtual bool update(const Vector & /* s  */, const Vector &  /* y */ ) = 0;


    virtual Matrix & get_Hessian() = 0;
    virtual Matrix & get_Hessian_inv() = 0;


private:
    Scalar num_tol_;
    bool initialized_;

};
}

#endif //UTOPIA_HESSIAN_APPROXIMATIONS_HPP
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
    
    // refresh vectors
    virtual bool update(const Vector & /* s  */, const Vector &  /* y */ ) = 0;    

    Scalar num_tol()const 
    {
        return num_tol_; 
    }

    void num_tol(Scalar & tol ) 
    {
        num_tol_ = tol; 
    }

    bool initialized() const 
    {
        return initialized_; 
    }
    
    std::function< void(const Vector &, Vector &) >  get_apply_Hinv()
    {
        std::function< void(const Vector &, Vector &) > my_func = 
        [this](const Vector &x, Vector & result)
            {
                this->apply_Hinv(x, result); 
            }; 

        return my_func; 
    }
    
    std::function< void(const Vector &, Vector &) >  get_compute_uHinvv_dot()
    {
        std::function< void(const Vector &, Vector &) > my_func = 
        [this](const Vector &x, Vector & result)
            {
                this->compute_uHinvv_dot(x, result); 
            }; 

        return my_func; 
    }


    std::function< void(const Vector &, Vector &) >  get_apply_H()
    {
        std::function< void(const Vector &, Vector &) > my_func = 
        [this](const Vector &x, Vector & result)
            {
                this->apply_H(x, result); 
            }; 

        return my_func; 
    }

    std::function< Scalar(const Vector &, const Vector &) >  get_compute_uHv_dot()
    {
        std::function< Scalar(const Vector &, const Vector &) > my_func = 
        [this](const Vector &x, const Vector & result)
            {
                return this->compute_uHv_dot(x, result); 
            }; 

        return my_func; 
    }


    std::function< Scalar(const Vector &) >  get_compute_uHu_dot()
    {
        std::function< Scalar(const Vector &) > my_func = 
        [this](const Vector &x)
            {
                return this->compute_uHu_dot(x); 
            }; 

        return my_func; 
    }    


protected: 
    void initialized(const bool init) 
    {
        initialized_ = init; 
    }

    // to be factored out 
    virtual bool constrained_solve(const Vector &/*x*/, const Vector &/*g*/, const Vector & /*lb*/, const Vector & /*ub*/, Vector & /*s*/, const Scalar & /*a  */) const {return false; }

    // applications of inverse of Hessian
    virtual bool apply_Hinv(const Vector & /* g */, Vector & /*s */) const  = 0;
    virtual Scalar compute_uHinvv_dot(const Vector &/*u*/, const Vector & /*v*/) const {return false; }

    // applications of Hessian
    virtual bool apply_H(const Vector & /*v*/ , Vector & /*r */) const  = 0 ; 
    virtual Scalar compute_uHv_dot(const Vector &/*u*/, const Vector & /*v*/) const {return 0; }
    virtual Scalar compute_uHu_dot(const Vector &/*u*/) const {return 0; }

    virtual Matrix & get_Hessian() = 0; 
    virtual Matrix & get_Hessian_inv() = 0; 


private:
    Scalar num_tol_;
    bool initialized_;

};
}

#endif //UTOPIA_HESSIAN_APPROXIMATIONS_HPP
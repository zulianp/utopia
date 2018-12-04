#ifndef UTOPIA_HESSIAN_APPROXIMATIONS_HPP
#define UTOPIA_HESSIAN_APPROXIMATIONS_HPP

#include "utopia_Core.hpp"
#include "utopia_TrivialPreconditioners.hpp"


namespace utopia
{

template<class Vector>
class HessianApproximation : public virtual Clonable 
{
    typedef UTOPIA_SCALAR(Vector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


public:

    HessianApproximation(): num_tol_(1e-12), initialized_(false)
    {
        
    }

    virtual ~HessianApproximation() { }

    virtual void initialize() = 0;
    virtual bool update(const Vector & /* s  */, const Vector &  /* y */ ) = 0;    
    virtual void reset() = 0; 

    virtual HessianApproximation<Vector> * clone() const override = 0;

    // applications of inverse of Hessian
    virtual bool apply_Hinv(const Vector & /* g */, Vector & /*s */) const  = 0;
    virtual bool apply_H(const Vector & /*v*/ , Vector & /*r */) const  = 0; 
    

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

    void initialized(const bool init) 
    {
        initialized_ = init; 
    }

    
    virtual Scalar compute_uHinvv_dot(const Vector & u, const Vector & v) const 
    {        
        Vector help; 
        this->apply_Hinv(v, help); 
        return dot(u, help); 
    }
    
    virtual Scalar compute_uHv_dot(const Vector & u , const Vector & v) const 
    {
        Vector help; 
        this->apply_H(v, help); 
        return dot(u, help); 
    }

    virtual Scalar compute_uHu_dot(const Vector & u) const
    {
        Vector help; 
        this->apply_H(u, help); 
        return dot(u, help); 
    }


    std::shared_ptr< FunctionOperator<Vector> > build_apply_Hinv()
    {
        std::function< void(const Vector &, Vector &) > my_func = 
        [this](const Vector &x, Vector & result)
        {
            this->apply_Hinv(x, result); 
        }; 

        return std::make_shared<FunctionOperator<Vector> >(my_func); 
    }



    std::shared_ptr< FunctionPreconditioner<Vector> > build_Hinv_precond()
    {
        std::function< void(const Vector &, Vector &) > my_func = 
        [this](const Vector &x, Vector & result)
        {
            this->apply_Hinv(x, result); 
        }; 

        return std::make_shared<FunctionPreconditioner<Vector> >(my_func); 
    }


    std::shared_ptr< FunctionOperator<Vector> > build_compute_uHinvv_dot()
    {
        std::function< Scalar(const Vector &, const Vector &) > my_func = 
        [this](const Vector &x, const Vector & result)
        {
            return this->compute_uHinvv_dot(x, result); 
        }; 

        return std::make_shared<FunctionOperator<Vector> >(my_func); 
    }


    std::shared_ptr< FunctionOperator<Vector> > build_apply_H()
    {
        std::function< void(const Vector &, Vector &) > my_func = 
        [this](const Vector &x, Vector & result)
        {
            this->apply_H(x, result); 
        }; 

        return std::make_shared<FunctionOperator<Vector> >(my_func); 
    }


    std::shared_ptr< FunctionOperator<Vector> > build_compute_uHv_dot()
    {
        std::function< Scalar(const Vector &, const Vector &) > my_func = 
        [this](const Vector &x, const Vector & result)
        {
            return this->compute_uHv_dot(x, result); 
        }; 

        return std::make_shared<FunctionOperator<Vector> >(my_func); 
    }



    std::shared_ptr< FunctionOperator<Vector> > build_compute_uHu_dot()
    {
        std::function< Scalar(const Vector &) > my_func = 
        [this](const Vector &x)
        {
            return this->compute_uHu_dot(x); 
        }; 

        return std::make_shared<FunctionOperator<Vector> >(my_func); 
    }


private:
    Scalar num_tol_;
    bool initialized_;

};
}

#endif //UTOPIA_HESSIAN_APPROXIMATIONS_HPP
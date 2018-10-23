#ifndef UTOPIA_MATRIX_FREE_SOLVER_INTERFACE_HPP
#define UTOPIA_MATRIX_FREE_SOLVER_INTERFACE_HPP

#include "utopia_Core.hpp"


namespace utopia
{

template<class Matrix, class Vector>
class MatrixFreeSolverInterface
{
    typedef UTOPIA_SCALAR(Vector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


public:

    MatrixFreeSolverInterface(): initialized_(false)
    {

    }

    virtual ~MatrixFreeSolverInterface() { }

    virtual void initialize(const std::shared_ptr <HessianApproximation<Matrix, Vector> > & approx)
    {
        apply_Hinv_         = approx->get_apply_Hinv(); 
        compute_uHinvv_dot_ = approx->get_compute_uHinvv_dot(); 

        apply_H_            = approx->get_apply_H(); 
        compute_uHv_dot_    = approx->get_compute_uHv_dot(); 
        compute_uHu_dot_    = approx->get_compute_uHu_dot(); 

        initialized_ = true; 
    }

    virtual  void initialized(const bool init) 
    {
        initialized_ = init; 
    }

    virtual bool initialized() const 
    {
        return initialized_; 
    }

    virtual  void apply_H(const Vector & v, Vector & r) const 
    {
        apply_H_(v, r); 
    }

    virtual void apply_Hinv(const Vector & v, Vector & r) const 
    {
        apply_Hinv_(v, r); 
    }

    virtual Scalar compute_uHu_dot(const Vector & v) const 
    {
        return compute_uHu_dot_(v); 
    }


private:
    bool initialized_;

    std::function< void(const Vector &, Vector &) > apply_Hinv_; 
    std::function< Scalar(const Vector &, Vector &) > compute_uHinvv_dot_; 
    
    std::function< void(const Vector &, Vector &) > apply_H_; 
    std::function< Scalar(const Vector &, const Vector &) > compute_uHv_dot_;
    std::function< Scalar(const Vector &) > compute_uHu_dot_;


};


// template<class Matrix, class Vector>
// class  QuasiDirectSolver:   public MatrixFreeSolverInterface<Matrix, Vector>, 
//                             public DirectSolver<Matrix, Vector> 
// {

//     public: 
//         virtual QuasiDirectSolver * clone() const override
//         {
//             return new QuasiDirectSolver(*this);
//         }


//         virtual bool apply(const Vector &rhs, Vector &sol) override
//         {
//             if(this->initialized())
//                 this->apply_Hinv(rhs, sol); 
//             else
//                 utopia_error("QuasiDirectSolver not initialized! \n"); 

//             return false; 
//         }        

// };

//////////////////////////////////////////////////////////////////////////////////////////////////////////

    template<class Vector>
    class FunctionOperator final: public Operator<Vector> 
    {
        public:
            FunctionOperator(const std::function< void(const Vector &, Vector &) > operator_action)
            : operator_action_(operator_action)
            {}

            bool apply(const Vector &rhs, Vector &ret) const override
            {
                operator_action_(rhs, ret); 
                return true;
            }


        private:
            std::function< void(const Vector &, Vector &) > operator_action_; 
    };



    // TODO:: do checks if action is inverse ...
    template<class Vector>
    class QuasiLinearSolver: public MatrixFreeLinearSolver<Vector>
    {
        public:
            virtual ~QuasiLinearSolver() {}

            virtual bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol)
            {
                A.apply(rhs, sol); 
                return true; 
            }

    };




}

#endif //UTOPIA_MATRIX_FREE_SOLVER_INTERFACE_HPP
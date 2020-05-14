#ifndef UTOPIA_HESSIAN_APPROXIMATIONS_HPP
#define UTOPIA_HESSIAN_APPROXIMATIONS_HPP

#include "utopia_Core.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include <memory>

namespace utopia
{

template<class Vector>
class HessianApproximation : public virtual Clonable, public virtual Configurable
{
    typedef UTOPIA_SCALAR(Vector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

    using Communicator = typename Traits<Vector>::Communicator;

public:

    class FunctionOperator final: public Operator<Vector>
    {
        public:

            FunctionOperator(
                HessianApproximation &parent,
                const std::function< void(const Vector &, Vector &) > operator_action)
            : parent(parent), operator_action_(operator_action)
            {}

            bool apply(const Vector &rhs, Vector &ret) const override
            {
                operator_action_(rhs, ret);
                return true;
            }

            inline Communicator &comm() override
            {
                return parent.comm();
            }

            inline const Communicator &comm() const override
            {
                return parent.comm();
            }

            inline Size size() const override
            {
                return parent.size();
            }

            inline Size local_size() const override
            {
                return parent.local_size();
            }

        private:
            HessianApproximation &parent;
            std::function< void(const Vector &, Vector &) > operator_action_;
    };

    HessianApproximation(): num_tol_(1e-12), initialized_(false)
    {}

    ~HessianApproximation() override {}

    virtual void initialize(const Vector &x_k, const Vector & /* g */)
    {
        comm_ = std::shared_ptr<Communicator>(x_k.comm().clone());

        size_.set_dims(2);
        size_.set(0, x_k.size());
        size_.set(1, x_k.size());

        local_size_.set_dims(2);
        local_size_.set(0, x_k.local_size());
        local_size_.set(1, x_k.local_size());
    }

    virtual bool update(const Vector & /* s  */, const Vector &  /* y */ , const Vector &  /* x */ , const Vector &  /* g */ ) = 0;
    virtual void reset() = 0;

    HessianApproximation<Vector> *clone() const override = 0;

    // applications of inverse of Hessian
    virtual bool apply_Hinv(const Vector & /* g */, Vector & /*s */)  = 0;
    virtual bool apply_H(const Vector & /*v*/ , Vector & /*r */)  = 0;

    void read(Input &in) override { in.get("num_tol", num_tol_); }

    void print_usage(std::ostream &os) const override {
        this->print_param_usage(os, "num_tol", "double", "Numerical tolerance.", "1e-12");
    }

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

    virtual Scalar compute_uHinvv_dot(const Vector & u, const Vector & v) 
    {
        Vector help;
        this->apply_Hinv(v, help);
        return dot(u, help);
    }

    virtual Scalar compute_uHv_dot(const Vector & u , const Vector & v)
    {
        Vector help;
        this->apply_H(v, help);
        return dot(u, help);
    }

    virtual Scalar compute_uHu_dot(const Vector & u)
    {
        Vector help;
        this->apply_H(u, help);
        return dot(u, help);
    }


    std::shared_ptr<FunctionOperator> build_apply_Hinv()
    {
        std::function< void(const Vector &, Vector &) > my_func =
        [this](const Vector &x, Vector & result)
        {
            this->apply_Hinv(x, result);
        };

        return std::make_shared<FunctionOperator>(*this, my_func);
        // return std::make_shared<FunctionOperator>(*this, &HessianApproximation::apply_Hinv);
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

    std::shared_ptr<FunctionOperator> build_compute_uHinvv_dot()
    {
        std::function< Scalar(const Vector &, const Vector &) > my_func =
        [this](const Vector &x, const Vector & result)
        {
            return this->compute_uHinvv_dot(x, result);
        };

        return std::make_shared<FunctionOperator>(*this, my_func);
    }

    std::shared_ptr<FunctionOperator> build_apply_H()
    {
        std::function< void(const Vector &, Vector &) > my_func =
        [this](const Vector &x, Vector & result)
        {
            this->apply_H(x, result);
        };

        return std::make_shared<FunctionOperator>(*this, my_func);
    }

    std::shared_ptr<FunctionOperator> build_compute_uHv_dot()
    {
        std::function< Scalar(const Vector &, const Vector &) > my_func =
        [this](const Vector &x, const Vector & result)
        {
            return this->compute_uHv_dot(x, result);
        };

        return std::make_shared<FunctionOperator>(*this, my_func);
    }

    std::shared_ptr<FunctionOperator> build_compute_uHu_dot()
    {
        std::function< Scalar(const Vector &) > my_func =
        [this](const Vector &x)
        {
            return this->compute_uHu_dot(x);
        };

        return std::make_shared<FunctionOperator>(*this, my_func);
    }

    virtual Size size() const
    {
        return size_;
    }

    virtual Size local_size() const
    {
        return local_size_;
    }

    virtual Communicator &comm()
    {
        assert(comm_);
        return *comm_;
    }

    virtual const Communicator &comm() const
    {
        assert(comm_);
        return *comm_;
    }

private:
    std::shared_ptr<Communicator> comm_;
    Scalar num_tol_;
    bool initialized_;
    Size size_, local_size_;
};
}

#endif //UTOPIA_HESSIAN_APPROXIMATIONS_HPP
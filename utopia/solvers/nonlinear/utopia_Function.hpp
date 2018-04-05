#ifndef UTOPIA_SOLVER_FUNCTION_HPP
#define UTOPIA_SOLVER_FUNCTION_HPP

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"

namespace utopia 
{
    /**
     * @brief      Base class for Nonlinear Function. All application context needed by solver is usually provided inside of this functions.
     *             In optimization settings, user needs to supply value(energy), gradient, hessian.
     *
     * @todo       Intorduce approximate Hessian updates strategies, e.g. BFGS, ... 
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class Function 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix); 

        virtual ~Function() { }

        virtual bool value(const Vector &/*point*/, Scalar &/*value*/) const = 0;
        virtual bool gradient(const Vector &/*point*/, Vector &/*result*/) const = 0;
        virtual bool hessian(const Vector &x, Matrix &H) const = 0;

        virtual bool hessian(const Vector &/*point*/, Matrix &/*result*/, Matrix &/*preconditioner*/) const 
        {
            return false;
        }

        virtual bool has_preconditioner() const {
            return false;
        }

        virtual bool update(const Vector &/*point*/) { return true; }
        virtual bool initialize_hessian(Matrix &/*H*/, Matrix &H_pre) const
        {
            return false;
        }

        /**
         * @brief Allows to solvers to reuse allocated vectors and matrices
         */
        class Data {
        public:
            std::shared_ptr<Matrix> H;
            std::shared_ptr<Matrix> H_pre;
            std::shared_ptr<Vector> g;

            void init()
            {
                if(!H) { H = std::make_shared<Matrix>(); }
                if(!H_pre) { H_pre = std::make_shared<Matrix>(); }
                if(!g) { g = std::make_shared<Vector>(); }
            }
        };
        
        inline std::shared_ptr<Data> data() const {
            return data_;
        }

        Function()
        : data_(std::make_shared<Data>())
        {}

    private:
        std::shared_ptr<Data> data_;

    };
}
#endif //UTOPIA_SOLVER_FUNCTION_HPP

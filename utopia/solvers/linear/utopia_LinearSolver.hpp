#ifndef UTOPIA_SOLVER_SOLVER_H
#define UTOPIA_SOLVER_SOLVER_H

#include <string>
#include "utopia_Core.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Utils.hpp"
#include "utopia_Clonable.hpp"

namespace  utopia
{
    /**
     * @brief      The base class for linear solvers.
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class LinearSolver : public Preconditioner<Vector> {
    public:
        typedef UTOPIA_SCALAR(Vector)           Scalar;

        virtual ~LinearSolver() {}

        virtual bool apply(const Vector &rhs, Vector &sol) override = 0;


        virtual void read(Input &in) override
        {
            Preconditioner<Vector>::read(in); 
        }
        
        virtual void print_usage(std::ostream &os) const override
        { 
            Preconditioner<Vector>::print_usage(os); 
        }

        /**
         * @brief      Solve routine.
         * @param[in]  A
         * @param[in]  b
         * @param      x0
         *
         * @return
         */
        virtual bool solve(const Matrix &A, const Vector &b, Vector &x0)
        {
            update(make_ref(A));
            return apply(b, x0);
        }


        /*! @brief if overriden the subclass has to also call this one first
         */
        virtual void update(const std::shared_ptr<const Matrix> &op)
        {
            op_ = op;
        }

        inline const std::shared_ptr<const Matrix> &get_operator() const
        {
            assert(op_);
            return op_;
        }

        inline bool has_operator() const
        {
            return static_cast<bool>(op_);
        }

        virtual LinearSolver * clone() const override = 0;
    private:
        std::shared_ptr<const Matrix> op_;
    };

}

#endif //UTOPIA_SOLVER_SOLVER_H

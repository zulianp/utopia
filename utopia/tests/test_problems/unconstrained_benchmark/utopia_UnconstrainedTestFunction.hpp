#ifndef UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS
#define UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"


namespace utopia
{

    // This function is used for implementation of test functions for unconstrained nonlinear benchmark
    template<class Matrix, class Vector>
    class UnconstrainedTestFunction : public Function<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        virtual ~UnconstrainedTestFunction() { }

        virtual Vector initial_guess() const = 0;
        virtual const Vector & exact_sol() const = 0;
        virtual Scalar min_function_value() const = 0;

        virtual std::string name() const = 0;
        virtual SizeType dim() const = 0;

        virtual bool exact_sol_known() const
        {
            return true;
        }

        virtual bool parallel() const
        {
            return false;
        }
    };
}
#endif //UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS

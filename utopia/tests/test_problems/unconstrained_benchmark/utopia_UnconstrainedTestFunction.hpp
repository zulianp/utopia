#ifndef UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS
#define UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"


namespace utopia
{

    // This function is used for implementation of test functions for unconstrained nonlinear benchmark
    template<class Matrix, class Vector>
    class UnconstrainedTestFunctionInterface
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        virtual ~UnconstrainedTestFunctionInterface() { }

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

        virtual void describe() const 
        {
            if(mpi_world_rank() == 0){
                std::cout<< name() <<",   Globalsize:" << dim() << ",   parallel:  " << parallel() <<  ",   sol. known:"<< exact_sol_known() << "  \n"; 
            }
        }

    };


    template<class Matrix, class Vector>
    class UnconstrainedTestFunction : public UnconstrainedTestFunctionInterface<Matrix, Vector>, public Function<Matrix, Vector>
    {
        public:
            DEF_UTOPIA_SCALAR(Matrix)
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        virtual ~UnconstrainedTestFunction() { }
    };

    template<class Matrix, class Vector>
    class UnconstrainedExtendedTestFunction : public UnconstrainedTestFunctionInterface<Matrix, Vector>, public ExtendedFunction<Matrix, Vector>
    {
        public:
            DEF_UTOPIA_SCALAR(Matrix)
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        virtual ~UnconstrainedExtendedTestFunction() { }
    };  

}
#endif //UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS

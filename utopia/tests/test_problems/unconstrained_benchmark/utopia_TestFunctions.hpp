#ifndef UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS
#define UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include <iomanip>
#include <locale>


namespace utopia
{

    // This function is used for implementation of test functions for unconstrained nonlinear benchmark
    template<class Matrix, class Vector>
    class TestFunctionInterface
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        virtual ~TestFunctionInterface() { }

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
                std::cout<< name() <<",   Globalsize:";  
                
                std::string numWithCommas = std::to_string(dim());
                int insertPosition = numWithCommas.length() - 3;
                while (insertPosition > 0) {
                    numWithCommas.insert(insertPosition, ",");
                    insertPosition-=3;
                }

                std::cout<<numWithCommas << ",   parallel:  " << parallel() <<  ",   sol. known:"<< exact_sol_known() << "  \n"; 
            }
        }
    };


    template<class Matrix, class Vector>
    class UnconstrainedTestFunction : virtual public TestFunctionInterface<Matrix, Vector>, virtual public Function<Matrix, Vector>
    {
        public:
            DEF_UTOPIA_SCALAR(Matrix);
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        virtual ~UnconstrainedTestFunction() { }
    };

    template<class Matrix, class Vector>
    class UnconstrainedExtendedTestFunction : virtual public TestFunctionInterface<Matrix, Vector>, virtual public ExtendedFunction<Matrix, Vector>
    {
        public:
            DEF_UTOPIA_SCALAR(Matrix);
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        virtual ~UnconstrainedExtendedTestFunction() { }
    };  

    template<class Matrix, class Vector>
    class ConstrainedTestFunction : virtual public TestFunctionInterface<Matrix, Vector>, virtual public Function<Matrix, Vector>
    {
        public:
            DEF_UTOPIA_SCALAR(Matrix);
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
            typedef utopia::BoxConstraints<Vector>                  BoxConstraints;

        virtual ~ConstrainedTestFunction() { }

        virtual const BoxConstraints & box_constraints() const
        {
            return constraints_;
        }

        virtual void set_box_constraints(const BoxConstraints & box)
        {
            constraints_ = box;
        }        

        virtual const Vector & upper_bound() const
        {
          if(!constraints_.upper_bound()){
            utopia_error("ConstrainedTestFunction::upper bound does not exist. \n");
          }

          return *constraints_.upper_bound();
        }

        virtual const Vector & lower_bound() const
        {
          if(!constraints_.lower_bound()){
            utopia_error("ConstrainedTestFunction::lower bound does not exist. \n");
          }

          return *constraints_.lower_bound();
        }

        virtual bool has_lower_bound() const
        {
          return constraints_.has_lower_bound();
        }

        virtual bool has_upper_bound() const
        {
          return constraints_.has_upper_bound();
        }


        protected:
            BoxConstraints                  constraints_;    

    };

    template<class Matrix, class Vector>
    class ConstrainedExtendedTestFunction : virtual public TestFunctionInterface<Matrix, Vector>, virtual public ExtendedFunction<Matrix, Vector>
    {
        public:
            DEF_UTOPIA_SCALAR(Matrix);
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        virtual ~ConstrainedExtendedTestFunction() { }

        virtual bool upper_bound(Vector &/*ub*/) const = 0;
        virtual bool lower_bound(Vector &/*lb*/) const = 0;   

        virtual bool has_upper_bound() const = 0;
        virtual bool has_lower_bound() const = 0;                     
    };      

}
#endif //UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS

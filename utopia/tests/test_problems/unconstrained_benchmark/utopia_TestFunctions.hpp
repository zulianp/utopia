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
        using Traits   = utopia::Traits<Vector>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

        virtual ~TestFunctionInterface() { }

        virtual Vector initial_guess() const= 0;
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

        virtual Vector initial_guess() const override= 0;
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

        virtual bool is_feasible(Vector & x)
        {
            if(!constraints_.has_upper_bound() && !constraints_.has_lower_bound())
                return true;

            if(empty(help_) || size(help_)!=size(x))
            {
                help_ = 0.0*x;
            }

            if(constraints_.has_upper_bound() && constraints_.has_lower_bound())
            {
                const auto &ub = *constraints_.upper_bound();
                const auto &lb = *constraints_.lower_bound();

                {
                    auto d_lb   = const_device_view(lb);
                    auto d_ub   = const_device_view(ub);
                    auto d_x    = const_device_view(x);

                    parallel_each_write(help_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        Scalar li = d_lb.get(i);
                        Scalar ui = d_ub.get(i);
                        Scalar xi = d_x.get(i);

                        return (xi < li || xi > ui) ? 1.0: 0.0;
                    });
                }

                return (sum(help_) > 0.0)? false : true;
            }
            else if(constraints_.has_upper_bound() && !constraints_.has_lower_bound())
            {
                const auto &ub = *constraints_.upper_bound();

                {
                    auto d_ub   = const_device_view(ub);
                    auto d_x    = const_device_view(x);

                    parallel_each_write(help_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        Scalar ui = d_ub.get(i);
                        Scalar xi = d_x.get(i);

                        return (xi > ui)? 1.0: 0.0;
                    });
                }

                return (sum(help_) > 0.0)? false : true;
            }
            else
            {
                const auto &lb = *constraints_.lower_bound();

                {
                    auto d_lb   = const_device_view(lb);
                    auto d_x    = const_device_view(x);

                    parallel_each_write(help_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        Scalar li = d_lb.get(i);
                        Scalar xi = d_x.get(i);

                        return (xi < li)? 1.0: 0.0;
                    });
                }

                return (sum(help_) > 0.0)? false : true;
            }
        }



        protected:
            BoxConstraints  constraints_;
            Vector          help_;

    };

    template<class Matrix, class Vector>
    class ConstrainedExtendedTestFunction : virtual public TestFunctionInterface<Matrix, Vector>, virtual public ExtendedFunction<Matrix, Vector>
    {
        public:
            DEF_UTOPIA_SCALAR(Matrix);
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
            typedef utopia::BoxConstraints<Vector>                  BoxConstraints;

        virtual ~ConstrainedExtendedTestFunction() { }

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

        virtual Vector initial_guess() const override= 0;

        virtual bool is_feasible(Vector & x)
        {
            if(!constraints_.has_upper_bound() && !constraints_.has_lower_bound())
                return true;

            if(empty(help_) || size(help_)!=size(x))
            {
                help_ = 0.0*x;
            }

            if(constraints_.has_upper_bound() && constraints_.has_lower_bound())
            {
                const auto &ub = *constraints_.upper_bound();
                const auto &lb = *constraints_.lower_bound();

                {
                    auto d_lb   = const_device_view(lb);
                    auto d_ub   = const_device_view(ub);
                    auto d_x    = const_device_view(x);

                    parallel_each_write(help_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        Scalar li = d_lb.get(i);
                        Scalar ui = d_ub.get(i);
                        Scalar xi = d_x.get(i);

                        return (xi < li || xi > ui) ? 1.0: 0.0;
                    });
                }

                return (sum(help_) > 0.0)? false : true;
            }
            else if(constraints_.has_upper_bound() && !constraints_.has_lower_bound())
            {
                const auto &ub = *constraints_.upper_bound();

                {
                    auto d_ub   = const_device_view(ub);
                    auto d_x    = const_device_view(x);

                    parallel_each_write(help_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        Scalar ui = d_ub.get(i);
                        Scalar xi = d_x.get(i);

                        return (xi > ui)? 1.0: 0.0;
                    });
                }

                return (sum(help_) > 0.0)? false : true;
            }
            else
            {
                const auto &lb = *constraints_.lower_bound();

                {
                    auto d_lb   = const_device_view(lb);
                    auto d_x    = const_device_view(x);

                    parallel_each_write(help_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        Scalar li = d_lb.get(i);
                        Scalar xi = d_x.get(i);

                        return (xi < li)? 1.0: 0.0;
                    });
                }

                return (sum(help_) > 0.0)? false : true;
            }
        }

        protected:
            BoxConstraints  constraints_;
            Vector          help_;
    };

}
#endif //UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS

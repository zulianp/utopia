#ifndef UTOPIA_QP_WITH_CONSTRAINTS
#define UTOPIA_QP_WITH_CONSTRAINTS

#include <vector>
#include <assert.h>
#include "utopia_TestFunctions.hpp"


namespace utopia
{

    /**
     * @brief      QP example with constraints
     */
    template<class Matrix, class Vector>
    class QPwithConstraints : public ConstrainedTestFunction<Matrix, Vector>
    {
    public:
        typedef UTOPIA_SIZE_TYPE(DVectord) SizeType;
        typedef UTOPIA_SCALAR(DVectord) Scalar;

        ~QPwithConstraints(){}


        QPwithConstraints()
        {
            const std::string data_path = Utopia::instance().get("data_path");
            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs_);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A_);

            Matrix I_1, I_2, I_3;

            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

             // from coarse to fine
            interpolation_operators_.push_back(std::move(I_1));
            interpolation_operators_.push_back(std::move(I_2));
            interpolation_operators_.push_back(std::move(I_3));


            lb_ = 0 * rhs_;                                           // lower bound
            ub_ = 0.1 * local_values(local_size(rhs_).get(0), 1);     // upper bound

        }

        bool value(const Vector &x, typename Vector::Scalar &energy) const override
        {
            energy = 0.5 * dot(x, A_ * x) - dot(x, rhs_);
            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            g = A_ * x - rhs_;
            return true;
        }

        bool hessian(const Vector & /*x*/, Matrix &H) const override
        {
            H = A_;
            return true;
        }

        bool upper_bound(Vector & upper, const Scalar & s) const
        {
            if(s==-999999)
                upper = ub_;
            else
                upper = s * local_values(local_size(rhs_).get(0), 1);
            return true;
        }

        bool lower_bound(Vector & lower, const Scalar & s) const
        {
            if(s==-999999)
                lower = lb_;
            else
                lower = s * local_values(local_size(rhs_).get(0), 1);
            return true;
        }

        bool upper_bound(Vector & upper) const override
        {
            upper = ub_;
            return true;
        }

        bool lower_bound(Vector & lower) const override
        {
            lower = lb_;
            return true;
        }

        bool has_upper_bound() const override
        {
            return true;
        }

        bool has_lower_bound() const override
        {
            return true;
        }


        std::vector<Matrix> & interpolation_operators()
        {
            return interpolation_operators_;
        }

        Vector & rhs()
        {
            return rhs_;
        }

        virtual Vector initial_guess() const override
        {   
            return (0 * rhs_);
        }
        
        virtual const Vector & exact_sol() const override
        {
            Vector empty; 
            return empty; 
        }
        

        virtual Scalar min_function_value() const override
        {   
            // not known
            return 0.0; 
        }

        virtual std::string name() const override
        {
            return "Laplace2D";
        }
        
        virtual SizeType dim() const override
        {
            return size(rhs_).get(0); 
        }

        virtual bool exact_sol_known() const override
        {
            return false;
        }

        virtual bool parallel() const override
        {
            return true;
        }


    private:
        Matrix A_;
        Vector rhs_;
        Vector ub_;
        Vector lb_;
        std::vector <Matrix> interpolation_operators_;


    };
}

#endif //UTOPIA_QP_WITH_CONSTRAINTS

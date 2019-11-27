#ifndef UTOPIA_POISSON_1D_HPP
#define UTOPIA_POISSON_1D_HPP

#include "utopia.hpp"
#include "utopia_TestFunctions.hpp"


namespace utopia
{

    template<typename Matrix, typename Vector>
    class Poisson1D final:  virtual public UnconstrainedExtendedTestFunction<Matrix, Vector>, 
                            virtual public ConstrainedExtendedTestFunction<Matrix, Vector>
    {
        public:
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
            typedef UTOPIA_SCALAR(Vector) Scalar;


        Poisson1D(const SizeType & n, const SizeType & problem_type=1):     pi_(3.14159265358979323846), 
                                                                            problem_type_(problem_type),  
                                                                            n_(n)
        { 

            if(problem_type_==1){
                assembly_problem_type1(); 
            }else if(problem_type_==2){
                assembly_problem_type2(); 
            }else if(problem_type_==3){
                assembly_problem_type3(); 
            }            
            else{
                utopia_error("Poisson1D:: problem type non-existent"); 
            }

            // this->constraints_ = make_box_constaints(std::make_shared<Vector>(values(n_, -9e9)),
            //                                          std::make_shared<Vector>(values(n_, 9e9)));    
        }


        bool initialize_hessian(Matrix &H, Matrix &/*H_pre*/) const override
        {
            H = H_;
            set_zero_rows(H, bc_indices_, 1.);

            return true;
        }

        bool value(const Vector &x, Scalar &value) const override
        {
            *A_help_  = (H_) * x; 
            value = 0.5 * dot(x, *A_help_) - dot(rhs_, x); 
            return true;
        }

        bool gradient_no_rhs(const Vector &x, Vector &g) const override
        {   
            g = (H_ * x) - rhs_; 

            {
                Write<Vector>   w(g); 
                Read<Vector>    read(x); 
                Read<Vector>    re(exact_sol_); 

                Range r = range(g);

                if(r.begin() == 0)  {
                    g.set(0, (exact_sol_.get(0) - x.get(0))); 
                }

                if(r.end() == n_)  {
                    g.set(n_-1, (exact_sol_.get(n_-1) - x.get(n_-1))); 
                }
            }

            return true;
        }        

        bool hessian(const Vector &x, Matrix &H) const override
        {
            H = H_; 
            set_zero_rows(H, bc_indices_, 1.);
            return true;
        }

        bool hessian(const Vector &x, Matrix &result, Matrix &prec) const override
        {
            UTOPIA_UNUSED(x);
            UTOPIA_UNUSED(result);
            UTOPIA_UNUSED(prec);
            return false;
        }

        bool has_preconditioner() const override
        {
            return false;
        }

        Vector initial_guess() const override
        {   
            return x0_; 
        }
        
        const Vector & exact_sol() const override
        {
            return exact_sol_; 
        }
        
        Scalar min_function_value() const override
        {   
            // depends on the solution to which we converged to 
            std::cout<<"Poisson1D:: min_function_value :: wrong.... \n"; 
            return -1.012; 
        }

        std::string name() const override
        {
            return "Poisson1D";
        }
        
        SizeType dim() const override
        {
            return n_; 
        }

        bool exact_sol_known() const override
        {
            return true;
        }

        bool parallel() const override
        {
            return true;
        }   


    private: 
        void assemble_laplacian_1D(Matrix &M)
        {
            {
                // n x n matrix with maximum 3 entries x row
                Write<Matrix> w(M);
                Range r = row_range(M);
                auto n = size(M).get(0);                

                for(SizeType i = r.begin(); i != r.end(); ++i) {
                    if(i > 0) {
                        M.set(i, i - 1, -1.0);
                    }

                    if(i < n-1) {
                        M.set(i, i + 1, -1.0);
                    }

                    if(i == 0 || i == n - 1) {
                        M.set(i, i, 1.);
                    } else {
                        M.set(i, i, 2.0);
                    }
                }
            }

            auto n = size(M).get(0);
            M *= 1./(h_*h_);  
        }

        void init_memory()
        {
            H_ = sparse(n_, n_, 3); 
            assemble_laplacian_1D(H_);            

            rhs_ = values(n_, 0.0); 
            x0_ = values(n_, 0.0); 
            exact_sol_ = values(n_, 0.0);
            A_help_ = make_unique<Vector>(values(n_, 0.0));
        }


        void assembly_problem_type1()
        {
            a_ = 0.0; 
            b_ = 2.0 * pi_;  

            L_ = b_ - a_; 
            h_ = L_ / (n_-1); 

            init_memory(); 

            {
                parallel_each_write(rhs_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    Scalar xi = (h_*i); 
                    if(i==0){
                        return xi * std::cos(xi); 
                    }
                    else if(i==n_-1){
                        return xi * std::cos(xi); 
                    }
                    else
                    {
                        // return (2.0* device::sin(xi)) + (xi*device::cos(xi)); 
                        return (2.0* std::sin(xi)) + (xi* std::cos(xi)); 
                    }
                });

                parallel_each_write(exact_sol_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    Scalar xi = (h_*i); 
                    // return xi * device::cos(xi); 
                    return xi * std::cos(xi); 
                });   

                parallel_each_write(x0_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    Scalar xi = (h_*i); 
                    if(i==0){
                        return xi * std::cos(xi); 
                    }
                    else if(i==n_-1){
                        return xi * std::cos(xi); 
                    }
                    else{
                        return 0.0; 
                    }
                });                                
            }


            Vector bc_markers = values(n_, 0.0);
            {
                Write<Vector> wv(bc_markers); 
                Range r = range(bc_markers);

                if(r.begin() == 0)  {
                    bc_markers.set(0, 1.0);
                }

                if(r.end() == n_)  {
                    bc_markers.set(n_-1, 1.0);
                }
            }

            bc_indices_.push_back(0.0);
            bc_indices_.push_back(n_-1);

            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, x0_);
        }

        void assembly_problem_type2()
        {
            a_ = 0.0; 
            b_ = 1.0;  

            L_ = b_ - a_; 
            h_ = L_ / (n_-1); 

            init_memory(); 

            {
                parallel_each_write(rhs_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    Scalar xi = (h_*i); 
                    if(i==0){
                        return 0.0;
                    }
                    else if(i==n_-1){
                        return 0.0;
                    }
                    else
                    {
                        return 10.0; 
                    }
                });
            }


            Vector bc_markers = values(n_, 0.0);
            {
                Write<Vector> wv(bc_markers); 
                Range r = range(bc_markers);

                if(r.begin() == 0)  {
                    bc_markers.set(0, 1.0);
                }

                if(r.end() == n_)  {
                    bc_markers.set(n_-1, 1.0);
                }
            }

            bc_indices_.push_back(0.0);
            bc_indices_.push_back(n_-1);

            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, x0_);
            ExtendedFunction<Matrix, Vector>::set_rhs(rhs_);


            Vector upper_bound = values(n_, 0.0); 
            {
                parallel_each_write(upper_bound, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    Scalar xi = (h_*i); 
                    return 0.5 + ((xi - 0.5) * (xi - 0.5));
                });                
            }            

            this->constraints_ = make_upper_bound_constraints(std::make_shared<Vector>(upper_bound)); 

        }        



        void assembly_problem_type3()
        {
            a_ = 0.0; 
            b_ = 1.0;  

            L_ = b_ - a_; 
            h_ = L_ / (n_-1); 

            init_memory(); 

            {
                parallel_each_write(rhs_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    Scalar xi = (h_*i); 
                    if(i==0){
                        return 0.0;
                    }
                    else if(i==n_-1){
                        return 0.0;
                    }
                    else
                    {
                        return 60.0; 
                    }
                });
            }


            Vector bc_markers = values(n_, 0.0);
            {
                Write<Vector> wv(bc_markers); 
                Range r = range(bc_markers);

                if(r.begin() == 0)  {
                    bc_markers.set(0, 1.0);
                }

                if(r.end() == n_)  {
                    bc_markers.set(n_-1, 1.0);
                }
            }

            bc_indices_.push_back(0.0);
            bc_indices_.push_back(n_-1);

            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, x0_);
            ExtendedFunction<Matrix, Vector>::set_rhs(rhs_);


            Vector upper_bound = values(n_, 0.0); 
            {
                parallel_each_write(upper_bound, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    Scalar xi = (h_*i); 
                    Scalar periods = 4.0; 
                    Scalar c = 2.0 * 3.14 * periods * (2.0*xi - 1.0)/2.0;
                    // device::cos
                    return  0.5 + ((std::cos(c)) - 0.5) * ((std::cos(c)) - 0.5);
                });                
            }            

            this->constraints_ = make_upper_bound_constraints(std::make_shared<Vector>(upper_bound)); 

        }   


    private: 
        const Scalar pi_; 
        const SizeType problem_type_; 

        Scalar a_, b_; 
        Scalar n_, L_, h_;         

        std::vector<SizeType> bc_indices_; 

        Matrix H_; 
        Vector rhs_;
        Vector x0_; 
        Vector exact_sol_; 

        std::unique_ptr<Vector>  A_help_; 


    }; 

}
#endif
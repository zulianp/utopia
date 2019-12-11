#ifndef UTOPIA_MOREBV_1D_HPP
#define UTOPIA_MOREBV_1D_HPP

#include "utopia.hpp"
#include "utopia_TestFunctions.hpp"


namespace utopia
{

    template<typename Matrix, typename Vector>
    class Morebv1D final:   virtual public UnconstrainedExtendedTestFunction<Matrix, Vector>, 
                            virtual public ConstrainedExtendedTestFunction<Matrix, Vector>
    {
        public:
            typedef UTOPIA_SIZE_TYPE(Vector)    SizeType;
            typedef UTOPIA_SCALAR(Vector)       Scalar;

        Morebv1D(const SizeType & n): n_(n)
        { 
            assembly_problem_type(); 
        }

        bool initialize_hessian(Matrix &H, Matrix &/*H_pre*/) const override
        {
            H = H_;
            set_zero_rows(H, bc_indices_, 1.);

            return true;
        }

        bool value(const Vector &x, Scalar &value) const override
        {
            A_help1_->set(0.0); 
            auto n = size(x).get(0);   

            {
                auto d_x  = const_device_view(x);
                auto d_t = const_device_view(coords_);

                parallel_each_write( *A_help1_ , UTOPIA_LAMBDA(const SizeType &i) -> Scalar 
                {
                    
                    Scalar xi = d_x.get(i);
                    Scalar ti = d_t.get(i);
                    Scalar element, item; 

                    if(i==0){
                        Scalar xi_p = d_x.get(i+1);
                        element = xi + h_ + 1.0;
                        item = 2. * xi - xi_p + 0.5*(h_*h_*std::pow(element, 3));
                        return item*item; 
                    }


                    if(i==n-1){
                        Scalar xi_m = d_x.get(i-1);
                        element = xi + ti + 1.0;
                        item = 2. * xi - xi_m + 0.5*(h_*h_*std::pow(element, 3));
                        return item*item; 
                    }

                    Scalar xi_m = d_x.get(i-1);
                    Scalar xi_p = d_x.get(i+1);
                    element = xi + ti + 1.0;
                    item = 2. * xi - xi_m - xi_p + 0.5*(h_*h_*std::pow(element, 3));

                    return item*item; 
                });
            }      

            value = sum(*A_help1_); 

            return true;
        }

        bool gradient_no_rhs(const Vector &x, Vector &g) const override
        {   
            if(empty(g)){
                g = local_zeros(local_size(x));
            }
            else{
                g.set(0.0);
            }

            auto n = size(x).get(0);   
            {
                Read<Vector>    d_x(x); 
                Read<Vector>    d_t(coords_); 
                Write<Vector>   w_t(g); 

                Range r = range(g);

                for(SizeType i = r.begin(); i != r.end(); ++i) 
                {
                    Scalar xi = x.get(i);
                    Scalar ti = coords_.get(i);
                    Scalar element, item; 

                    if(i==0){
                        Scalar xi_p = x.get(i+1);
                        element = xi + h_ + 1.0;
                        item = 2. * xi - xi_p + 0.5*(h_*h_*std::pow(element, 3));
                        g.add(i,    2.0*item *(2.0+ h_*h_ *1.5*std::pow(element,2))); 
                        g.add(i+1, -2.0*item); 
                    }
                    else if(i==n-1){
                        Scalar xi_m = x.get(i-1);
                        element = xi + ti + 1.0;
                        item = 2. * xi - xi_m + 0.5*(h_*h_*std::pow(element, 3));
                        g.add(i,    4.0*item *(2.0+ h_*h_ *1.5*std::pow(element,2))); 
                        g.add(i-1, -2.0*item); 
                    }
                    else
                    {
                        Scalar xi_m = x.get(i-1);
                        Scalar xi_p = x.get(i+1);
                        element = xi + ti + 1.0;
                        item = 2. * xi - xi_m - xi_p + 0.5*(h_*h_*std::pow(element, 3));

                        g.add(i-1,  -2.0*item);
                        g.add(i,    2.0*item*(2.0+ h_*h_ * 1.5 * std::pow(element, 2)));
                        g.add(i+1,  -2.0*item);
                    }
                }

                if(r.begin() == 0){                
                    g.set(0,  -x.get(0)); 
                }
                if(r.end() == n_){                
                    g.set(n_-1, -x.get(n_-1)); 
                }                
            }      

            return true;
        }    


        bool hessian(const Vector &x, Matrix &H) const override
        {
            if(empty(H)){
                H = H_; 
            }

            H = local_identity(local_size(H_).get(0), local_size(H_).get(1)); 


            // // (d)
            // *A_help1_ = (H_) * x; 
            // {
            //     Write<Matrix> w(H);
            //     Read<Vector> r1(*A_help1_);
            //     Read<Vector> r2(coords_);
            //     Read<Vector> r3(x);

            //     Range r = row_range(H);
            //     auto n = size(H).get(0);                

            //     Scalar hh = h_*h_; 

            //     for(SizeType i = r.begin(); i != r.end(); ++i) 
            //     {
            //         if(i > 0){
            //             // Scalar tx1  = coords_.get(i) + x.get(i) + 1.0; 
            //             // Scalar eb   =  (3./2. * hh * tx1*tx1 + 2.0);                        
            //             Scalar v = -1.0; // * 1.0; 
            //             H.set(i, i - 1, v);
            //         }

            //         if(i < n-1){
            //             // Scalar tx1  = coords_.get(i) + x.get(i) + 1.0; 
            //             // Scalar eb   =  (3./2. * hh * tx1*tx1 + 2.0);                        
            //             Scalar v = -1.0; // * 1.0;                         
            //             H.set(i, i + 1, v);
            //         }

            //         if(i == 0 || i == n - 1){
            //             H.set(i, i, 1.0);
            //         }
            //         else{
            //             Scalar tx1  = coords_.get(i) + x.get(i) + 1.0; 
            //             Scalar eb   =  (3./2. * hh * tx1*tx1 + 2.0);                        
            //             Scalar v = 0.5 * hh * tx1*tx1*tx1; 
            //             v = (v + A_help1_->get(i)) * 6.0 * hh * tx1; 
            //             v = v + (2.0*eb*eb); 

            //             H.set(i, i, v);
            //         }
            //     }
            // }


            // set_zero_rows(H, bc_indices_, 1.);

            // disp(H); 
            // exit(0); 



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
            return 0.0; 
        }

        std::string name() const override
        {
            return "Morebv1D";
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
        }

        void init_memory()
        {
            H_ = sparse(n_, n_, 3); 
            assemble_laplacian_1D(H_);            

            coords_     = values(n_, 0.0); 
            ones_       = values(n_, 1.0); 
            twos_       = values(n_, 2.0); 
            x0_         = values(n_, 0.0); 
            exact_sol_  = values(n_, 0.0);
            A_help1_    = make_unique<Vector>(values(n_, 0.0));
            A_help2_    = make_unique<Vector>(values(n_, 0.0));
        }


        void assembly_problem_type()
        {
            a_ = 0.0; 
            b_ = 1.0;  

            L_ = b_ - a_; 
            h_ = L_ / (n_-1); 

            init_memory(); 

            {
                parallel_each_write(coords_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    return (h_*i);
                });

                // see More, Garbbow, Hillstrom
                parallel_each_write(x0_, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    Scalar xi = (h_*i); 
                    if(i==0){
                        return 0.0;
                    }
                    else if(i==n_-1){
                        return 0.0;
                    }
                    else{
                        return xi * (xi - 1.0); 
                    }
                });     

            }

            Vector bc_markers = values(n_, 0.0);
            {
                Write<Vector> wv(bc_markers); 
                Range r = range(bc_markers);

                if(r.begin() == 0)  {
                    bc_markers.set(0, 1.0);
                    bc_indices_.push_back(0.0);
                }

                if(r.end() == n_)  {
                    bc_markers.set(n_-1, 1.0);
                    bc_indices_.push_back(n_-1);
                }
            }


            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, x0_);
            

            // Vector upper_bound = values(n_, 0.0); 
            // {
            //     parallel_each_write(upper_bound, UTOPIA_LAMBDA(const SizeType i) -> Scalar
            //     {
            //         Scalar xi = (h_*i); 
            //         return 0.5 + ((xi - 0.5) * (xi - 0.5));
            //     });                
            // }            

            // this->constraints_ = make_upper_bound_constraints(std::make_shared<Vector>(upper_bound)); 

        }        


    private: 
        Scalar a_, b_; 
        Scalar n_, L_, h_;         

        std::vector<SizeType> bc_indices_; 

        Matrix H_; 
        Vector coords_;
        Vector ones_; 
        Vector twos_; 
        Vector x0_; 
        Vector exact_sol_; 


        std::unique_ptr<Vector>  A_help1_; 
        std::unique_ptr<Vector>  A_help2_; 

    }; 

}
#endif